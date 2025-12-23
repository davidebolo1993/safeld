#include "memory_pool.h"
#include "utils.h"

VectorPool g_vector_pool(2000);

VectorPool::VectorPool(size_t max_size) : max_size_(max_size) {
    // Silent initialization - no logging here
}

// NO DESTRUCTOR - let compiler generate default one

std::vector<double> VectorPool::acquireDoubleVector(size_t reserve_size) {
    std::lock_guard<std::mutex> lock(double_mutex_);
    
    if (!double_vectors_.empty()) {
        auto vec = std::move(double_vectors_.front());
        double_vectors_.pop();
        vec.clear();
        if (reserve_size > 0) {
            vec.reserve(reserve_size);
        }
        return vec;
    }
    
    std::vector<double> vec;
    if (reserve_size > 0) {
        vec.reserve(reserve_size);
    }
    return vec;
}

std::vector<std::string> VectorPool::acquireStringVector(size_t reserve_size) {
    std::lock_guard<std::mutex> lock(string_mutex_);
    
    if (!string_vectors_.empty()) {
        auto vec = std::move(string_vectors_.front());
        string_vectors_.pop();
        vec.clear();
        if (reserve_size > 0) {
            vec.reserve(reserve_size);
        }
        return vec;
    }
    
    std::vector<std::string> vec;
    if (reserve_size > 0) {
        vec.reserve(reserve_size);
    }
    return vec;
}

void VectorPool::releaseDoubleVector(std::vector<double>&& vec) {
    std::lock_guard<std::mutex> lock(double_mutex_);
    
    if (double_vectors_.size() < max_size_) {
        double_vectors_.push(std::move(vec));
    }
}

void VectorPool::releaseStringVector(std::vector<std::string>&& vec) {
    std::lock_guard<std::mutex> lock(string_mutex_);
    
    if (string_vectors_.size() < max_size_) {
        string_vectors_.push(std::move(vec));
    }
}

size_t VectorPool::getDoubleVectorPoolSize() const {
    std::lock_guard<std::mutex> lock(double_mutex_);
    return double_vectors_.size();
}

size_t VectorPool::getStringVectorPoolSize() const {
    std::lock_guard<std::mutex> lock(string_mutex_);
    return string_vectors_.size();
}

PooledDoubleVector acquireDoubleVector(size_t reserve_size) {
    return PooledDoubleVector(g_vector_pool.acquireDoubleVector(reserve_size), &g_vector_pool);
}

std::vector<std::string> acquireStringVector(size_t reserve_size) {
    return g_vector_pool.acquireStringVector(reserve_size);
}

void releaseStringVector(std::vector<std::string>&& vec) {
    g_vector_pool.releaseStringVector(std::move(vec));
}

