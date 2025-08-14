#pragma once

#include <vector>
#include <memory>
#include <random>
#include <future>
#include <queue>
#include <mutex>
#include <condition_variable>
#include "vcf_processor.h"

struct ProcessedVariant {
    size_t index;
    std::string id;
    std::string chrom;
    int pos;
    std::string ref;
    std::string alt;
    std::vector<double> synthetic_dosages;
};

class SimulationEngine {
private:
    int n_traits_;
    int n_samples_;
    int n_workers_;
    std::vector<std::vector<double>> traits_matrix_;
    
    // Thread pool components
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex queue_mutex_;
    std::condition_variable cv_;
    bool stop_;
    
    void initializeTraitsMatrix();
    void workerThread();
    ProcessedVariant processVariant(const VariantRecord& variant, size_t index);
    
public:
    SimulationEngine(int n_traits, int n_samples, int n_workers = 0);
    ~SimulationEngine();
    
    void initialize();
    std::vector<ProcessedVariant> simulateVariants(
        const std::vector<std::unique_ptr<VariantRecord>>& variants);
    
    // Getters
    int getNumTraits() const { return n_traits_; }
    int getNumSamples() const { return n_samples_; }
};

