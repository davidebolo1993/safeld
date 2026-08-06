#pragma once

#include <string>
#include <vector>
#include <memory>
#include <chrono>

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

void setVerboseLogging(bool enabled);
bool isVerboseLogging();
void logInfo(const std::string& message);
void logDebug(const std::string& message);
void logWarning(const std::string& message);
void logError(const std::string& message);
