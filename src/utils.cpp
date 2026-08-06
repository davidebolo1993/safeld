#include "utils.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <mutex>
#include <atomic>

static std::mutex log_mutex;
static std::atomic<bool> verbose_logging{false};

void setVerboseLogging(bool enabled) {
    verbose_logging.store(enabled);
}

bool isVerboseLogging() {
    return verbose_logging.load();
}

void logInfo(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cout << "[INFO] " << message << std::endl;
}

void logDebug(const std::string& message) {
    if (!verbose_logging.load()) {
        return;
    }
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cout << "[DEBUG] " << message << std::endl;
}

void logWarning(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cout << "[WARN] " << message << std::endl;
}

void logError(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cerr << "[ERROR] " << message << std::endl;
}

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
    double elapsed_time = elapsed();
    logDebug(name + " completed in " + std::to_string(elapsed_time) + " seconds");
}

void Timer::reset() {
    start_time = std::chrono::high_resolution_clock::now();
}

double Timer::elapsed() const {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    return duration.count() / 1000000.0;
}
