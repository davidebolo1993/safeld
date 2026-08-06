#include "utils.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <mutex>
#include <atomic>
#include <ctime>
#include <unistd.h>

namespace {

std::mutex log_mutex;
std::atomic<bool> verbose_logging{false};

// Guarded by log_mutex; stages are sequential so contention is irrelevant.
std::string log_module = "safeld";

// True when a progress bar currently owns the last line of the terminal, so the
// next log line can wipe it before printing.
bool bar_on_screen = false;

std::string timestamp() {
    std::time_t now = std::time(nullptr);
    std::tm tm_buf{};
    localtime_r(&now, &tm_buf);
    char buf[16];
    std::strftime(buf, sizeof(buf), "%H:%M:%S", &tm_buf);
    return std::string(buf);
}

// Erase a progress bar that is sitting on the current line. Must be called with
// log_mutex held.
void clearBarLocked() {
    if (bar_on_screen) {
        std::cerr << "\r\033[K" << std::flush;
        bar_on_screen = false;
    }
}

void emit(std::ostream& os, const char* level, const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    clearBarLocked();
    os << "[" << timestamp() << "] [" << log_module << "] "
       << level << " " << message << std::endl;
}

bool stderrIsTerminal() {
    return isatty(fileno(stderr)) != 0;
}

std::string renderBarBody(double fraction, int width) {
    if (fraction < 0.0) fraction = 0.0;
    if (fraction > 1.0) fraction = 1.0;
    int filled = static_cast<int>(fraction * width + 0.5);
    std::string body;
    body.reserve(width * 3);
    for (int i = 0; i < width; ++i) {
        body += (i < filled) ? "█" : "░";
    }
    return body;
}

}  // namespace

void setVerboseLogging(bool enabled) {
    verbose_logging.store(enabled);
}

bool isVerboseLogging() {
    return verbose_logging.load();
}

void setLogModule(const std::string& module) {
    std::lock_guard<std::mutex> lock(log_mutex);
    log_module = module;
}

std::string currentLogModule() {
    std::lock_guard<std::mutex> lock(log_mutex);
    return log_module;
}

LogModule::LogModule(const std::string& module) {
    previous_ = currentLogModule();
    setLogModule(module);
}

LogModule::~LogModule() {
    setLogModule(previous_);
}

void logInfo(const std::string& message) {
    emit(std::cout, "INFO ", message);
}

void logDebug(const std::string& message) {
    if (!verbose_logging.load()) {
        return;
    }
    emit(std::cout, "DEBUG", message);
}

void logWarning(const std::string& message) {
    emit(std::cout, "WARN ", message);
}

void logError(const std::string& message) {
    emit(std::cerr, "ERROR", message);
}

std::string formatCount(long long value) {
    bool negative = value < 0;
    unsigned long long v = negative ? static_cast<unsigned long long>(-value)
                                    : static_cast<unsigned long long>(value);
    std::string digits = std::to_string(v);
    std::string out;
    int count = 0;
    for (auto it = digits.rbegin(); it != digits.rend(); ++it) {
        if (count && count % 3 == 0) out.push_back(',');
        out.push_back(*it);
        count++;
    }
    if (negative) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

std::string formatBytes(unsigned long long bytes) {
    static const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB"};
    double value = static_cast<double>(bytes);
    int unit = 0;
    while (value >= 1024.0 && unit < 5) {
        value /= 1024.0;
        unit++;
    }
    std::ostringstream os;
    os << std::fixed << std::setprecision(value < 10.0 && unit > 0 ? 1 : 0)
       << value << " " << units[unit];
    return os.str();
}

std::string formatDuration(double seconds) {
    std::ostringstream os;
    if (seconds < 60.0) {
        os << std::fixed << std::setprecision(1) << seconds << "s";
    } else if (seconds < 3600.0) {
        int m = static_cast<int>(seconds) / 60;
        int s = static_cast<int>(seconds) % 60;
        os << m << "m" << std::setw(2) << std::setfill('0') << s << "s";
    } else {
        int h = static_cast<int>(seconds) / 3600;
        int m = (static_cast<int>(seconds) % 3600) / 60;
        os << h << "h" << std::setw(2) << std::setfill('0') << m << "m";
    }
    return os.str();
}

// ---------------------------------------------------------------------------
// ProgressBar
// ---------------------------------------------------------------------------

ProgressBar::ProgressBar(const std::string& label, long long total)
    : label_(label), total_(total), current_(0),
      interactive_(stderrIsTerminal()), finished_(false),
      start_(std::chrono::steady_clock::now()), last_draw_(start_) {
    if (total_ <= 0) {
        // Nothing to track; behave as a no-op so callers need no special case.
        finished_ = true;
    }
}

ProgressBar::~ProgressBar() {
    if (!finished_) {
        finish();
    }
}

void ProgressBar::update(long long current) {
    if (finished_) return;
    current_ = current;
    render(false);
}

void ProgressBar::increment(long long delta) {
    if (finished_) return;
    current_ += delta;
    render(false);
}

void ProgressBar::render(bool final_frame) {
    auto now = std::chrono::steady_clock::now();

    if (!final_frame) {
        // Redraw at most ~10x a second on a terminal; in a log file only on
        // decile boundaries, so a redirected run stays readable.
        if (interactive_) {
            if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_draw_).count() < 100) {
                return;
            }
        } else {
            // Redirected output: a heartbeat every 30s. Short stages then print
            // nothing but their completion line, and a long one still shows it
            // is alive without burying the log.
            if (std::chrono::duration_cast<std::chrono::seconds>(now - last_draw_).count() < 30) {
                return;
            }
        }
    }
    last_draw_ = now;

    double fraction = total_ > 0 ? static_cast<double>(current_) / static_cast<double>(total_) : 1.0;
    if (fraction > 1.0) fraction = 1.0;
    if (fraction < 0.0) fraction = 0.0;
    double elapsed = std::chrono::duration<double>(now - start_).count();

    if (interactive_) {
        std::ostringstream os;
        os << "\r\033[K  " << label_ << " " << renderBarBody(fraction, 24) << " "
           << std::fixed << std::setprecision(0) << (fraction * 100.0) << "% "
           << formatCount(current_) << "/" << formatCount(total_);
        if (fraction > 0.02 && !final_frame) {
            double eta = elapsed / fraction - elapsed;
            os << "  eta " << formatDuration(eta);
        }
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << os.str() << std::flush;
        bar_on_screen = true;
    } else if (!final_frame) {
        std::ostringstream os;
        os << label_ << ": " << std::fixed << std::setprecision(0) << (fraction * 100.0)
           << "% (" << formatCount(current_) << "/" << formatCount(total_) << ")";
        logInfo(os.str());
    }
}

void ProgressBar::finish(const std::string& suffix) {
    if (finished_) return;
    finished_ = true;

    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_).count();

    {
        std::lock_guard<std::mutex> lock(log_mutex);
        clearBarLocked();
    }

    std::ostringstream os;
    os << label_ << ": " << formatCount(current_) << "/" << formatCount(total_)
       << " in " << formatDuration(elapsed);
    if (!suffix.empty()) os << " (" << suffix << ")";
    logInfo(os.str());
}

// ---------------------------------------------------------------------------
// ProgressCounter
// ---------------------------------------------------------------------------

ProgressCounter::ProgressCounter(const std::string& label, const std::string& unit,
                                 long long report_every)
    : label_(label), unit_(unit), report_every_(report_every > 0 ? report_every : 100000),
      current_(0), last_reported_(0), interactive_(stderrIsTerminal()), finished_(false),
      start_(std::chrono::steady_clock::now()), last_draw_(start_) {
}

void ProgressCounter::increment(long long delta) {
    if (finished_) return;
    current_ += delta;

    if (interactive_) {
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_draw_).count() < 100) {
            return;
        }
        last_draw_ = now;
        double elapsed = std::chrono::duration<double>(now - start_).count();
        std::ostringstream os;
        os << "\r\033[K  " << label_ << " " << formatCount(current_) << " " << unit_;
        if (elapsed > 0.5) {
            os << "  (" << formatCount(static_cast<long long>(current_ / elapsed)) << "/s)";
        }
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << os.str() << std::flush;
        bar_on_screen = true;
    } else if (current_ - last_reported_ >= report_every_) {
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now - last_draw_).count() >= 30) {
            last_draw_ = now;
            last_reported_ = current_;
            logInfo(label_ + ": " + formatCount(current_) + " " + unit_);
        }
    }
}

void ProgressCounter::finish(const std::string& suffix) {
    if (finished_) return;
    finished_ = true;

    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_).count();

    {
        std::lock_guard<std::mutex> lock(log_mutex);
        clearBarLocked();
    }

    std::ostringstream os;
    os << label_ << ": " << formatCount(current_) << " " << unit_
       << " in " << formatDuration(elapsed);
    if (!suffix.empty()) os << " (" << suffix << ")";
    logInfo(os.str());
}

// ---------------------------------------------------------------------------

std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

double parseDouble(const std::string& str, double defaultValue) {
    try {
        return std::stod(str);
    } catch (const std::exception&) {
        return defaultValue;
    }
}

bool standardize(const std::vector<double>& data, std::vector<double>& out) {
    if (data.empty()) {
        return false;
    }

    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();

    // Mean-imputed entries sit exactly at the mean and so contribute nothing to
    // the sum of squares while still counting towards the denominator. This is
    // the usual convention (plink does the same) and shrinks sigma by
    // sqrt(n_observed / n) at variants with missing calls.
    double sq_sum = 0.0;
    for (double value : data) {
        sq_sum += (value - mean) * (value - mean);
    }
    double std_dev = std::sqrt(sq_sum / data.size());

    if (!(std_dev > 0.0)) {
        return false;
    }

    out.clear();
    out.reserve(data.size());

    for (double value : data) {
        out.push_back((value - mean) / std_dev);
    }

    return true;
}

std::vector<double> scaleToDosageRange(const std::vector<double>& data) {
    if (data.empty()) {
        return {};
    }

    auto minmax = std::minmax_element(data.begin(), data.end());
    double min_val = *minmax.first;
    double max_val = *minmax.second;

    double range = max_val - min_val;
    if (range == 0.0) {
        return std::vector<double>(data.size(), 1.0);
    }

    std::vector<double> scaled;
    scaled.reserve(data.size());

    for (double value : data) {
        double scaled_value = 2.0 * (value - min_val) / range;
        scaled.push_back(scaled_value);
    }

    return scaled;
}

Timer::Timer(const std::string& timer_name) : name(timer_name) {
    reset();
}

Timer::~Timer() {
    logDebug(name + " completed in " + formatDuration(elapsed()));
}

void Timer::reset() {
    start_time = std::chrono::high_resolution_clock::now();
}

double Timer::elapsed() const {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    return duration.count() / 1000000.0;
}
