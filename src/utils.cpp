#include "utils.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <mutex>

// Thread-safe logging
static std::mutex log_mutex;

void logInfo(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cout << "[INFO] " << message << std::endl;
}

void logWarning(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cout << "[WARN] " << message << std::endl;
}

void logError(const std::string& message) {
    std::lock_guard<std::mutex> lock(log_mutex);
    std::cerr << "[ERROR] " << message << std::endl;
}

// String utilities
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

// Statistical functions
std::vector<double> standardize(const std::vector<double>& data) {
    if (data.empty()) {
        return {};
    }
    
    // Calculate mean
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();
    
    // Calculate standard deviation
    double sq_sum = 0.0;
    for (double value : data) {
        sq_sum += (value - mean) * (value - mean);
    }
    double std_dev = std::sqrt(sq_sum / data.size());
    
    if (std_dev == 0.0) {
        return std::vector<double>(data.size(), 0.0);
    }
    
    // Standardize
    std::vector<double> standardized;
    standardized.reserve(data.size());
    
    for (double value : data) {
        standardized.push_back((value - mean) / std_dev);
    }
    
    return standardized;
}

std::vector<double> scaleToDosageRange(const std::vector<double>& data) {
    if (data.empty()) {
        return {};
    }
    
    // Find min and max
    auto minmax = std::minmax_element(data.begin(), data.end());
    double min_val = *minmax.first;
    double max_val = *minmax.second;
    
    double range = max_val - min_val;
    if (range == 0.0) {
        return std::vector<double>(data.size(), 1.0); // Default to middle of [0, 2]
    }
    
    // Scale to [0, 2] range
    std::vector<double> scaled;
    scaled.reserve(data.size());
    
    for (double value : data) {
        double scaled_value = 2.0 * (value - min_val) / range;
        scaled.push_back(scaled_value);
    }
    
    return scaled;
}

// Timer implementation
Timer::Timer(const std::string& timer_name) : name(timer_name) {
    reset();
}

Timer::~Timer() {
    double elapsed_time = elapsed();
    logInfo(name + " completed in " + std::to_string(elapsed_time) + " seconds");
}

void Timer::reset() {
    start_time = std::chrono::high_resolution_clock::now();
}

double Timer::elapsed() const {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    return duration.count() / 1000000.0; // Convert to seconds
}

