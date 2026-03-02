#pragma once

#include <vector>
#include <memory>
#include <mutex>
#include <queue>
#include <atomic>

template<typename T>
class ObjectPool {
private:
    std::queue<std::unique_ptr<T>> pool_;
    std::mutex mutex_;
    std::atomic<size_t> created_objects_{0};
    std::atomic<size_t> reused_objects_{0};
    size_t max_size_;
    
public:
    explicit ObjectPool(size_t max_size = 1000) : max_size_(max_size) {}
    
    template<typename... Args>
    std::unique_ptr<T> acquire(Args&&... args) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        if (!pool_.empty()) {
            auto obj = std::move(pool_.front());
            pool_.pop();
            reused_objects_++;
            return obj;
        }
        
        created_objects_++;
        return std::make_unique<T>(std::forward<Args>(args)...);
    }
    
    void release(std::unique_ptr<T> obj) {
        if (!obj) return;
        
        std::lock_guard<std::mutex> lock(mutex_);
        
        if (pool_.size() < max_size_) {
            pool_.push(std::move(obj));
        }
    }
    
    size_t getCreatedCount() const { return created_objects_.load(); }
    size_t getReusedCount() const { return reused_objects_.load(); }
    size_t getCurrentPoolSize() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return pool_.size();
    }
};

class VectorPool {
private:
    std::queue<std::vector<double>> double_vectors_;
    std::queue<std::vector<std::string>> string_vectors_;
    mutable std::mutex double_mutex_;
    mutable std::mutex string_mutex_;
    size_t max_size_;
    
public:
    explicit VectorPool(size_t max_size = 1000);
    
    std::vector<double> acquireDoubleVector(size_t reserve_size = 0);
    void releaseDoubleVector(std::vector<double>&& vec);
    
    std::vector<std::string> acquireStringVector(size_t reserve_size = 0);
    void releaseStringVector(std::vector<std::string>&& vec);
    
    size_t getDoubleVectorPoolSize() const;
    size_t getStringVectorPoolSize() const;
};

template<typename T>
class PooledObject {
private:
    std::unique_ptr<T> obj_;
    ObjectPool<T>* pool_;
    
public:
    PooledObject(std::unique_ptr<T> obj, ObjectPool<T>* pool)
        : obj_(std::move(obj)), pool_(pool) {}
    
    ~PooledObject() {
        if (obj_ && pool_) {
            pool_->release(std::move(obj_));
        }
    }
    
    PooledObject(const PooledObject&) = delete;
    PooledObject& operator=(const PooledObject&) = delete;
    
    PooledObject(PooledObject&& other) noexcept
        : obj_(std::move(other.obj_)), pool_(other.pool_) {
        other.pool_ = nullptr;
    }
    
    PooledObject& operator=(PooledObject&& other) noexcept {
        if (this != &other) {
            if (obj_ && pool_) {
                pool_->release(std::move(obj_));
            }
            obj_ = std::move(other.obj_);
            pool_ = other.pool_;
            other.pool_ = nullptr;
        }
        return *this;
    }
    
    T* operator->() { return obj_.get(); }
    const T* operator->() const { return obj_.get(); }
    T& operator*() { return *obj_; }
    const T& operator*() const { return *obj_; }
    T* get() { return obj_.get(); }
    const T* get() const { return obj_.get(); }
    
    bool valid() const { return obj_ != nullptr; }
};

class PooledDoubleVector {
private:
    std::vector<double> vec_;
    VectorPool* pool_;
    bool released_;
    
public:
    PooledDoubleVector(std::vector<double>&& vec, VectorPool* pool)
        : vec_(std::move(vec)), pool_(pool), released_(false) {}
    
    ~PooledDoubleVector() {
        if (!released_ && pool_) {
            pool_->releaseDoubleVector(std::move(vec_));
        }
    }
    
    PooledDoubleVector(const PooledDoubleVector&) = delete;
    PooledDoubleVector& operator=(const PooledDoubleVector&) = delete;
    
    PooledDoubleVector(PooledDoubleVector&& other) noexcept
        : vec_(std::move(other.vec_)), pool_(other.pool_), released_(other.released_) {
        other.released_ = true;
    }
    
    PooledDoubleVector& operator=(PooledDoubleVector&& other) noexcept {
        if (this != &other) {
            if (!released_ && pool_) {
                pool_->releaseDoubleVector(std::move(vec_));
            }
            vec_ = std::move(other.vec_);
            pool_ = other.pool_;
            released_ = other.released_;
            other.released_ = true;
        }
        return *this;
    }
    
    std::vector<double>& get() { return vec_; }
    const std::vector<double>& get() const { return vec_; }
    std::vector<double>& operator*() { return vec_; }
    const std::vector<double>& operator*() const { return vec_; }
    
    void resize(size_t size) { vec_.resize(size); }
    void reserve(size_t size) { vec_.reserve(size); }
    void clear() { vec_.clear(); }
    size_t size() const { return vec_.size(); }
    bool empty() const { return vec_.empty(); }
    
    double& operator[](size_t index) { return vec_[index]; }
    const double& operator[](size_t index) const { return vec_[index]; }
};

extern VectorPool g_vector_pool;
extern ObjectPool<std::vector<double>> g_double_vector_pool;
extern ObjectPool<std::vector<std::string>> g_string_vector_pool;

PooledDoubleVector acquireDoubleVector(size_t reserve_size = 0);
std::vector<std::string> acquireStringVector(size_t reserve_size = 0);
void releaseStringVector(std::vector<std::string>&& vec);
