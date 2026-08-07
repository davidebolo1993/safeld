#pragma once

#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <cstdint>

std::vector<std::string> split(const std::string& str, char delimiter);
double parseDouble(const std::string& str, double defaultValue = 0.0);

// Standardises `data` into `out`. Returns false, leaving `out` untouched, when
// the input is empty or has zero variance: such a vector cannot be standardised,
// and emitting a constant row instead would travel silently through the pipeline
// and surface as a synthetic dosage of exactly 1.0 for every trait.
bool standardize(const std::vector<double>& data, std::vector<double>& out);

std::vector<double> scaleToDosageRange(const std::vector<double>& data);

class Timer {
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::string name;

public:
    explicit Timer(const std::string& timer_name);
    ~Timer();
    void reset();
    double elapsed() const;
};

// ---------------------------------------------------------------------------
// Logging
//
// Every line is written as:
//     [HH:MM:SS] [module] LEVEL  message
// so that a log can be read back long after the run without guessing which
// stage produced which line. The module is a process-wide label set once per
// stage; see LogModule below.
// ---------------------------------------------------------------------------

void setVerboseLogging(bool enabled);
bool isVerboseLogging();

// Sets the label shown in the second bracket. Prefer the RAII LogModule.
void setLogModule(const std::string& module);
std::string currentLogModule();

// Scoped module label, restored on destruction.
class LogModule {
public:
    explicit LogModule(const std::string& module);
    ~LogModule();
private:
    std::string previous_;
};

void logInfo(const std::string& message);
void logDebug(const std::string& message);
void logWarning(const std::string& message);
void logError(const std::string& message);

// Human-readable helpers for counts and byte sizes, so messages stay readable
// at genomic scale ("1,234,567" rather than "1234567", "2.3 GB" rather than
// a byte count).
std::string formatCount(long long value);
std::string formatBytes(unsigned long long bytes);
std::string formatDuration(double seconds);

// ---------------------------------------------------------------------------
// Progress reporting
//
// On a terminal this draws a single self-overwriting bar. When output is
// redirected to a file it degrades to occasional one-line updates instead, so
// logs do not fill with control characters. Use it only where the total is
// known in advance; for open-ended scans use ProgressCounter.
// ---------------------------------------------------------------------------

class ProgressBar {
public:
    ProgressBar(const std::string& label, long long total);
    ~ProgressBar();

    void update(long long current);
    void increment(long long delta = 1);
    // Clears the bar and prints a completion line with the elapsed time.
    void finish(const std::string& suffix = "");

private:
    void render(bool final_frame);

    std::string label_;
    long long total_;
    long long current_;
    bool interactive_;
    bool finished_;
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point last_draw_;
};

// Open-ended counterpart for streams whose length is unknown until the end.
class ProgressCounter {
public:
    ProgressCounter(const std::string& label, const std::string& unit, long long report_every = 100000);
    void increment(long long delta = 1);
    void finish(const std::string& suffix = "");

private:
    std::string label_;
    std::string unit_;
    long long report_every_;
    long long current_;
    long long last_reported_;
    bool interactive_;
    bool finished_;
    std::chrono::steady_clock::time_point start_;
    std::chrono::steady_clock::time_point last_draw_;
};
