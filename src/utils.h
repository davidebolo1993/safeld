#pragma once

#include <string>
#include <vector>
#include <memory>
#include <chrono>

// Utility functions
std::vector<std::string> split(const std::string& str, char delimiter);
double parseDouble(const std::string& str, double defaultValue = 0.0);
std::vector<double> standardize(const std::vector<double>& data);
std::vector<double> scaleToDosageRange(const std::vector<double>& data);

// Timer class for performance monitoring
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

// Thread-safe logging
void logInfo(const std::string& message);
void logWarning(const std::string& message);
void logError(const std::string& message);

