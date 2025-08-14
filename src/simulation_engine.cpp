#include "simulation_engine.h"
#include "utils.h"
#include <algorithm>
#include <numeric>
#include <thread>
#include <cblas.h>

SimulationEngine::SimulationEngine(int n_traits, int n_samples, int n_workers)
    : n_traits_(n_traits), n_samples_(n_samples), stop_(false) {
    
    if (n_workers <= 0) {
        n_workers_ = std::thread::hardware_concurrency();
        if (n_workers_ == 0) n_workers_ = 4;  // Fallback
    } else {
        n_workers_ = n_workers;
    }
    
    logInfo("Simulation engine initialized with " + std::to_string(n_traits_) + 
            " traits, " + std::to_string(n_samples_) + " samples, " +
            std::to_string(n_workers_) + " workers");
}

SimulationEngine::~SimulationEngine() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        stop_ = true;
    }
    cv_.notify_all();
    
    for (auto& worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

void SimulationEngine::initializeTraitsMatrix() {
    Timer timer("Traits matrix generation");
    
    traits_matrix_.resize(n_traits_);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    
    for (int i = 0; i < n_traits_; ++i) {
        traits_matrix_[i].resize(n_samples_);
        for (int j = 0; j < n_samples_; ++j) {
            traits_matrix_[i][j] = dist(gen);
        }
    }
    
    logInfo("Generated synthetic traits matrix: " + std::to_string(n_traits_) + 
            " x " + std::to_string(n_samples_));
}

void SimulationEngine::workerThread() {
    while (true) {
        std::function<void()> task;
        
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            cv_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
            
            if (stop_ && tasks_.empty()) {
                return;
            }
            
            task = std::move(tasks_.front());
            tasks_.pop();
        }
        
        task();
    }
}

void SimulationEngine::initialize() {
    Timer timer("Simulation engine initialization");
    
    initializeTraitsMatrix();
    
    // Start worker threads
    workers_.reserve(n_workers_);
    for (int i = 0; i < n_workers_; ++i) {
        workers_.emplace_back(&SimulationEngine::workerThread, this);
    }
    
    logInfo("Started " + std::to_string(n_workers_) + " worker threads");
}

ProcessedVariant SimulationEngine::processVariant(const VariantRecord& variant, size_t index) {
    // Standardize dosages
    auto standardized = standardize(variant.dosages);
    
    // Compute synthetic dosages using matrix multiplication
    std::vector<double> synthetic_dosages(n_traits_);
    
    for (int t = 0; t < n_traits_; ++t) {
        double sum = cblas_ddot(n_samples_, standardized.data(), 1, 
                               traits_matrix_[t].data(), 1);
        synthetic_dosages[t] = sum / n_samples_;
    }
    
    // Scale to dosage range [0, 2]
    auto scaled = scaleToDosageRange(synthetic_dosages);
    
    ProcessedVariant processed;
    processed.index = index;
    processed.id = variant.id;
    processed.chrom = variant.chrom;
    processed.pos = variant.pos;
    processed.ref = variant.ref;
    processed.alt = variant.alt;
    processed.synthetic_dosages = std::move(scaled);
    
    return processed;
}

std::vector<ProcessedVariant> SimulationEngine::simulateVariants(
    const std::vector<std::unique_ptr<VariantRecord>>& variants) {
    
    Timer timer("Variant simulation");
    
    logInfo("Starting simulation of " + std::to_string(variants.size()) + " variants");
    
    std::vector<ProcessedVariant> results;
    results.reserve(variants.size());
    
    // Process variants sequentially to avoid thread exhaustion
    for (size_t i = 0; i < variants.size(); ++i) {
        auto processed = processVariant(*variants[i], i);
        results.push_back(std::move(processed));
        
        if (i % 1000 == 0 && i > 0) {
            logInfo("Simulated " + std::to_string(i) + "/" + 
                    std::to_string(variants.size()) + " variants");
        }
    }
    
    logInfo("Simulation completed for " + std::to_string(results.size()) + " variants");
    return results;
}


