#include "memory_pool.h"
#include "utils.h"

// Global pool instances
VectorPool g_vector_pool(2000);
ObjectPool<std::vector<double>> g_double_vector_pool(1000);
ObjectPool<std::vector<std::string>> g_string_vector_pool(500);

VectorPool::VectorPool(size_t max_size) : max_size_(max_size) {
    logInfo("Vector pool initialized with max size: " + std::to_string(max_size));
}

std::vector<double> VectorPool::acquireDoubleVector(size_t reserve_size) {
    std::lock_guard<std::mutex> lock(double_mutex_);
    
    if (!double_vectors_.empty()) {
        auto vec = std::move(double_vectors_.front());
        double_vectors_.pop();
        
        vec.clear();
        if (reserve_size > 0 && vec.capacity() < reserve_size) {
            vec.reserve(reserve_size);
        }
        
        return vec;
    }
    
    // Create new vector
    std::vector<double> vec;
    if (reserve_size > 0) {
        vec.reserve(reserve_size);
    }
    return vec;
}

void VectorPool::releaseDoubleVector(std::vector<double>&& vec) {
    std::lock_guard<std::mutex> lock(double_mutex_);
    
    if (double_vectors_.size() < max_size_) {
        vec.clear();  // Clear contents but keep capacity
        double_vectors_.push(std::move(vec));
    }
    // If pool is full, let the vector be destroyed
}

std::vector<std::string> VectorPool::acquireStringVector(size_t reserve_size) {
    std::lock_guard<std::mutex> lock(string_mutex_);
    
    if (!string_vectors_.empty()) {
        auto vec = std::move(string_vectors_.front());
        string_vectors_.pop();
        
        vec.clear();
        if (reserve_size > 0 && vec.capacity() < reserve_size) {
            vec.reserve(reserve_size);
        }
        
        return vec;
    }
    
    // Create new vector
    std::vector<std::string> vec;
    if (reserve_size > 0) {
        vec.reserve(reserve_size);
    }
    return vec;
}

void VectorPool::releaseStringVector(std::vector<std::string>&& vec) {
    std::lock_guard<std::mutex> lock(string_mutex_);
    
    if (string_vectors_.size() < max_size_) {
        vec.clear();  // Clear contents but keep capacity
        string_vectors_.push(std::move(vec));
    }
    // If pool is full, let the vector be destroyed
}

size_t VectorPool::getDoubleVectorPoolSize() const {
    std::lock_guard<std::mutex> lock(double_mutex_);
    return double_vectors_.size();
}

size_t VectorPool::getStringVectorPoolSize() const {
    std::lock_guard<std::mutex> lock(string_mutex_);
    return string_vectors_.size();
}

// Convenience functions
PooledDoubleVector acquireDoubleVector(size_t reserve_size) {
    auto vec = g_vector_pool.acquireDoubleVector(reserve_size);
    return PooledDoubleVector(std::move(vec), &g_vector_pool);
}

std::vector<std::string> acquireStringVector(size_t reserve_size) {
    return g_vector_pool.acquireStringVector(reserve_size);
}

void releaseStringVector(std::vector<std::string>&& vec) {
    g_vector_pool.releaseStringVector(std::move(vec));
}

