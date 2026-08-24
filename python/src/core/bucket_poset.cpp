
// MSVC: explicit standard includes (not pulled in transitively).
#include <string>
#include <string_view>
#include <memory>
#include <utility>
#include <cstdint>
// ================================================================
// File: bucketPOSet.cpp
// ================================================================

#include "bucket_poset.h"

#include <stdexcept>
#include <vector>

// ***********************************************
// Build
// ***********************************************

std::unique_ptr<BucketPOSet> BucketPOSet::Build(
                                                const std::vector<std::string>& elements,
                                                const std::vector<std::pair<std::string, std::string>>& comparabilities) {
    
    // Inizializzazione con ownership esclusiva
    std::unique_ptr<BucketPOSet> result(new BucketPOSet(elements.size()));
    
    // Richiama la funzione base ereditata da POSet
    result->FillBaseAttributes(elements, comparabilities, nullptr, true);
    result->BuildBucketOrder();
    
    // Semantica di spostamento automatica (RVO/NRVO)
    return result;
}

// ***********************************************
// BuildBucketOrder
// ***********************************************

void BucketPOSet::BuildBucketOrder() {
    const std::uint64_t n = this->size();
    
    for (std::uint64_t e_value = 0; e_value < n; ++e_value) {
        bool bucket_trovato = false;
        
        // HPC: Il ciclo itera sui Bucket
        for (auto it = buckets.begin(); it != buckets.end(); ++it) {
            auto bucket = *it;
            bool incomparabile = true;
            
            for (auto e_in_b : *bucket) {
                // Interroghiamo l'infrastruttura di base ereditata
                if (this->IsLessOrEqual(e_value, e_in_b) || this->IsLessOrEqual(e_in_b, e_value)) {
                    incomparabile = false;
                    break;
                }
            }
            
            if (incomparabile) {
                bucket->insert(e_value);
                bucket_trovato = true;
                break;
            }
        }
        
        if (!bucket_trovato) {
            auto new_bucket = std::make_shared<Bucket>(n);
            new_bucket->insert(e_value);
            buckets.push_back(std::move(new_bucket));
        }
    }
}

// ***********************************************
// Clone
// ***********************************************

std::shared_ptr<POSet> BucketPOSet::Clone() const {
    auto elements = std::vector<std::string>();
    elements.reserve(this->size());
    
    for (std::uint64_t i = 0; i < this->size(); ++i) {
        elements.emplace_back(this->GetElementName(i));
    }
    
    
    std::vector<std::pair<std::string, std::string>> comparabilities;
    comparabilities.reserve(this->size() * 4);
    
    for (std::uint64_t d1 = 0; d1 < this->size(); ++d1) {
        std::string_view v1_view = this->GetElementName(d1);
        
        for (std::uint64_t d2 = 0; d2 < this->size(); ++d2) {
            if (d1 != d2 && this->IsLessOrEqual(d1, d2)) {
                std::string_view v2_view = this->GetElementName(d2);
                comparabilities.emplace_back(std::string{v1_view}, std::string{v2_view});
            }
        }
    }
    
    return BucketPOSet::Build(elements, comparabilities);
}

// ***********************************************
// Getters and Replace
// ***********************************************

std::shared_ptr<Bucket> BucketPOSet::BucketAt(std::uint64_t layer) {
    return buckets.at(layer);
}

std::uint64_t BucketPOSet::LSize(std::uint64_t layer) const {
    return buckets.at(layer)->size();
}

std::uint64_t BucketPOSet::LSize() const {
    return buckets.size();
}

bool BucketPOSet::IsTotalOrder() const {
    return buckets.size() == this->size();
}

void BucketPOSet::ReplaceBuckets(BucketOrder& new_buckets) {
    this->buckets = std::move(new_buckets);
}

// ***********************************************
// Evaluation MRP
// ***********************************************

Tensor<double, 2> BucketPOSet::evaluationMRP() {
    return buckets.evaluationMRP();
}
