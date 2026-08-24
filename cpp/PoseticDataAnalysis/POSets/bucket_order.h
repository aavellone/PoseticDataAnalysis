// ================================================================
// File: bucketOrder.h
// ================================================================

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>
#include <cstddef>

#include "bucket.h"
#include "tensor.h"
#include "loss_function_mrp.h"

#include <vector>
#include <memory>
#include <string>
#include <tuple>
#include <cstdint>
#include <format>
#include <iterator>
#include <concepts>
#include <stdexcept>


// ----------------------------------------------------------------
// Alias Type HPC (shared_ptr per Matrice rimossi)
// ----------------------------------------------------------------

using BUCKETPTR = std::shared_ptr<Bucket>;

using BJRES = std::tuple<
    BUCKETPTR,
    Tensor<double, 2>,
    double
>;

using BSRES = std::tuple<
    BUCKETPTR,
    BUCKETPTR,
    Tensor<double, 2>,
    double
>;

class BucketOrder;

using RANKINGEL = std::tuple<
    std::shared_ptr<BucketOrder>,
    Tensor<double, 2>,
    double,
    std::string
>;

using RANKING = std::vector<RANKINGEL>;

// ----------------------------------------------------------------
// Class BucketOrder
// ----------------------------------------------------------------

class BucketOrder {
public:
    using bo_type        = std::vector<std::shared_ptr<Bucket>>;
    using iterator       = bo_type::iterator;
    using const_iterator = bo_type::const_iterator;
    
private:
    bo_type buckets_;
    std::uint64_t element_number_;
    
public:
    explicit BucketOrder(std::uint64_t e_number) : element_number_(e_number) {
        buckets_.reserve(32);
    }
    
    [[nodiscard]] std::uint64_t eSize() const noexcept { return element_number_; }
    [[nodiscard]] std::size_t size() const noexcept { return buckets_.size(); }
    [[nodiscard]] bool empty() const noexcept { return buckets_.empty(); }
    
    [[nodiscard]] iterator begin() { return buckets_.begin(); }
    [[nodiscard]] iterator end() { return buckets_.end(); }
    [[nodiscard]] const_iterator begin() const { return buckets_.begin(); }
    [[nodiscard]] const_iterator end() const { return buckets_.end(); }
    
    [[nodiscard]] iterator iter_at(std::uint64_t layer) {
        if (layer >= buckets_.size()) throw std::out_of_range("BucketOrder: layer out of range");
        auto it = buckets_.begin();
        std::advance(it, static_cast<std::ptrdiff_t>(layer));
        return it;
    }
    
    [[nodiscard]] const_iterator iter_at(std::uint64_t layer) const {
        if (layer >= buckets_.size()) throw std::out_of_range("BucketOrder: layer out of range");
        auto it = buckets_.begin();
        std::advance(it, static_cast<std::ptrdiff_t>(layer));
        return it;
    }
    
    // --- RISOLTO L'ERRORE CONST-CORRECTNESS QUI ---
    [[nodiscard]] std::shared_ptr<Bucket> at(std::uint64_t layer) { return *iter_at(layer); }
    [[nodiscard]] std::shared_ptr<Bucket> at(std::uint64_t layer) const { return *iter_at(layer); }
    
    [[nodiscard]] std::shared_ptr<Bucket> layer_at(std::uint64_t layer) { return *iter_at(layer); }
    [[nodiscard]] std::shared_ptr<Bucket> layer_at(std::uint64_t layer) const { return *iter_at(layer); }
    // ----------------------------------------------
    
    void push_back(std::shared_ptr<Bucket> b) { buckets_.push_back(std::move(b)); }
    void clear() { buckets_.clear(); }
    
    [[nodiscard]] std::shared_ptr<BucketOrder> copy() const {
        auto bo = std::make_shared<BucketOrder>(element_number_);
        for (const auto& b : buckets_) {
            auto new_b = std::make_shared<Bucket>(element_number_);
            for (auto v : *b) new_b->insert(v);
            bo->push_back(new_b);
        }
        return bo;
    }
    
    void addAt(const Bucket& elements_to_add, std::uint64_t layer) {
        auto b = this->at(layer);
        for (auto e : elements_to_add) b->insert(e);
    }
    
    void eraseAt(const Bucket& elements_to_remove, std::uint64_t layer) {
        auto b = this->at(layer);
        auto new_b = b->set_difference(elements_to_remove);
        *b = std::move(new_b);
    }
    
    void replace(const Bucket& new_bucket_data, std::uint64_t layer) {
        auto new_b = std::make_shared<Bucket>(element_number_);
        for(auto e : new_bucket_data) new_b->insert(e);
        *iter_at(layer) = new_b;
    }
    
    void insert(std::shared_ptr<Bucket> new_bucket, std::uint64_t layer) {
        buckets_.insert(iter_at(layer), std::move(new_bucket));
    }
    
    // ================================================================
    // Evaluation MRP
    // ================================================================
    
    [[nodiscard]] Tensor<double, 2> evaluationMRP() const {
        std::uint64_t totale_elementi = 0;
        for (auto bucket : buckets_) {
            totale_elementi += bucket->size();
        }
        
        Tensor<double, 2> mrp({totale_elementi, totale_elementi}, 0.0);
        std::vector<bool> used(totale_elementi, false);
        
        for (auto bucket : buckets_) {
            for (auto riga : *bucket) {
                mrp(riga, riga) = 1.0;
                
                for (std::uint64_t colonna = 0; colonna < totale_elementi; ++colonna) {
                    if (riga == colonna) continue;
                    
                    if (!used[colonna]) {
                        if (!bucket->contains(colonna)) {
                            mrp(riga, colonna) = 1.0;
                        }
                    }
                }
            }
            
            for (auto element : *bucket) {
                used[element] = true;
            }
        }
        
        return mrp;
    }
    
    // ================================================================
    // Metodi Algoritmici Template
    // ================================================================
    
    template <std::invocable<std::uint64_t, std::uint64_t> LeqFunc>
    static bool BSplitPreserve(const Bucket& b_down, const Bucket& b_up, LeqFunc is_leq) {
        for (auto v1 : b_up) {
            for (auto v2 : b_down) {
                if (is_leq(v1, v2)) return false;
            }
        }
        return true;
    }
    
    template <std::invocable<std::uint64_t, std::uint64_t> LeqFunc>
    static bool BJoinPreserve(const Bucket& b_down, const Bucket& b_up, const Bucket& b_top, LeqFunc is_leq) {
        for (auto v1 : b_up) {
            for (auto v2 : b_down) {
                if (is_leq(v1, v2)) return false;
            }
        }
        for (auto v1 : b_up) {
            for (auto v2 : b_top) {
                if (is_leq(v1, v2)) return false;
            }
        }
        return true;
    }
    
    template <std::invocable<std::uint64_t, std::uint64_t> LeqStart,
    std::invocable<std::uint64_t, std::uint64_t> LeqPres>
    BJRES BestJoinAt(LeqStart is_leq_start, LeqPres is_leq_pres, LossFunctionMRP& lfmrp, std::uint64_t bucket_position) {
        std::uint64_t bucket_order_size = this->eSize();
        double best_lfmrp_value = 0.0;
        Tensor<double, 2> best_mrp({bucket_order_size, bucket_order_size}, 0.0);
        BUCKETPTR best_b = nullptr;
        
        auto bucket = this->at(bucket_position);
        
        if (bucket->size() <= 1) {
            return std::make_tuple(nullptr, Tensor<double, 2>({bucket_order_size, bucket_order_size}, 0.0), 0.0);
        }
        
        auto b1 = std::make_shared<Bucket>(bucket_order_size);
        
        for (auto it_bucket = bucket->begin(); it_bucket != bucket->end() && b1->size() < bucket->size() - 1; ++it_bucket) {
            bool all_incomparable = true;
            auto bin = this->at(bucket_position);
            
            if (!all_incomparable) continue;
            
            b1->insert(*it_bucket);
            auto b2 = bucket->set_difference(*b1);
            
            if (!BucketOrder::BJoinPreserve(b2, *b1, *(this->at(bucket_position + 1)), is_leq_pres)) {
                continue;
            }
            
            auto bucket_order_copy = this->copy();
            bucket_order_copy->addAt(*b1, bucket_position + 1);
            bucket_order_copy->eraseAt(*b1, bucket_position);
            
            Tensor<double, 2> mrp = bucket_order_copy->evaluationMRP();
            auto lfmrp_value = lfmrp(mrp);
            
            if (lfmrp_value < best_lfmrp_value || best_b == nullptr) {
                best_lfmrp_value = lfmrp_value;
                best_mrp = std::move(mrp);
                best_b = std::make_shared<Bucket>(bucket_order_size);
                for(auto e : *b1) best_b->insert(e);
            }
        }
        
        return std::make_tuple(best_b, best_mrp, best_lfmrp_value);
    }
    
    template <std::invocable<std::uint64_t, std::uint64_t> LeqStart,
    std::invocable<std::uint64_t, std::uint64_t> LeqPres>
    BSRES BestSplitAt(LeqStart is_leq_start, LeqPres is_leq_pres, LossFunctionMRP& lfmrp, std::uint64_t bucket_position) {
        std::uint64_t bucket_order_size = this->eSize();
        double best_lfmrp_value = 0.0;
        Tensor<double, 2> best_mrp({bucket_order_size, bucket_order_size}, 0.0);
        BUCKETPTR best_b1 = nullptr;
        BUCKETPTR best_b2 = nullptr;
        
        auto bucket = this->at(bucket_position);
        if (bucket->size() <= 1) {
            return std::make_tuple(nullptr, nullptr, Tensor<double, 2>({bucket_order_size, bucket_order_size}, 0.0), 0.0);
        }
        
        auto b1 = std::make_shared<Bucket>(bucket_order_size);
        for (auto it_bucket = bucket->begin(); it_bucket != bucket->end() && b1->size() < bucket->size() - 1; ++it_bucket) {
            b1->insert(*it_bucket);
            auto b2 = bucket->set_difference(*b1);
            
            if (!BucketOrder::BSplitPreserve(*b1, b2, is_leq_pres)) {
                continue;
            }
            
            auto bucket_order_copy = this->copy();
            bucket_order_copy->replace(b2, bucket_position);
            bucket_order_copy->insert(b1, bucket_position);
            
            Tensor<double, 2> mrp = bucket_order_copy->evaluationMRP();
            auto lfmrp_value = lfmrp(mrp);
            
            if (lfmrp_value < best_lfmrp_value || best_b1 == nullptr) {
                best_lfmrp_value = lfmrp_value;
                best_mrp = std::move(mrp);
                
                best_b1 = std::make_shared<Bucket>(bucket_order_size);
                for(auto e : *b1) best_b1->insert(e);
                
                best_b2 = std::make_shared<Bucket>(bucket_order_size);
                for(auto e : b2) best_b2->insert(e);
            }
        }
        
        return std::make_tuple(best_b1, best_b2, best_mrp, best_lfmrp_value);
    }
};
