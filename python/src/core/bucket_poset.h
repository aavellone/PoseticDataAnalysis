// ================================================================
// File: bucketPOSet.h
// ================================================================

#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>

#include "poset.h"
#include "bucket_order.h"
#include "tensor.h"

#include <vector>
#include <list>
#include <string>
#include <memory>
#include <cstdint>



class BucketPOSet : public POSet {
private:
    BucketOrder buckets;
    
    explicit BucketPOSet(std::uint_fast64_t esize) : POSet(), buckets(esize) {}
    
public:
    [[nodiscard]] static std::unique_ptr<BucketPOSet> Build(
                    const std::vector<std::string>& elements,
                    const std::vector<std::pair<std::string, std::string>>& comparabilities);
    
    void ReplaceBuckets(BucketOrder& new_buckets);
    
    ~BucketPOSet() override = default;
    
    /**
     * @brief Calcola la matrice MRP del POSet interfacciandosi con BucketOrder.
     * Restituisce by-value per sfruttare RVO (no allocazioni dinamiche).
     */
    [[nodiscard]] Tensor<double, 2> evaluationMRP();
    
    bool IsTotalOrder() const override;
    std::uint_fast64_t LSize(std::uint_fast64_t layer) const;
    std::uint_fast64_t LSize() const;
    std::shared_ptr<Bucket> BucketAt(std::uint_fast64_t layer);
    
    BucketOrder& Buckets() {
        return this->buckets;
    }
    
    std::shared_ptr<POSet> Clone() const override;
    
private:
    void BuildBucketOrder();
};
