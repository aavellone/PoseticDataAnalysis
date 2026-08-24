#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>
#include <cstdint>
#include <functional>

#include <vector>
#include <string>
#include <numeric>
#include <limits>
#include <algorithm>
#include <ranges>
#include "linear_extension.h"
#include "bit_set.h"

using IdealSet = BitSet;

class TreeOfIdeals {
public:
    static constexpr std::uint64_t NO_LABEL = std::numeric_limits<std::uint64_t>::max();
    
private:
    std::vector<std::uint64_t> parent_;
    std::vector<std::uint64_t> labels_;
    std::vector<std::vector<std::uint64_t>> children_;
    std::vector<IdealSet> ideals_;
    
    // I predecessori ora sono nativamente immagazzinati in un array di BitSet
    std::vector<BitSet> im_pred_;
    LinearExtension le_for_conversion_;
    std::uint64_t root_;
    
public:
    // Il costruttore ora accetta direttamente std::vector<BitSet>
    TreeOfIdeals(std::vector<BitSet> im_pred,
                 LinearExtension le)
    : im_pred_(std::move(im_pred)), le_for_conversion_(std::move(le)) {
        
        std::uint64_t n = im_pred_.size();
        
        IdealSet ideal(n);
        ideal.FillUpTo(n);
        
        // Indicizzazione a base 0 corretta
        root_ = Left(static_cast<std::int64_t>(n) - 1, std::move(ideal));
    }
    
    // Getter rinominati secondo lo standard PascalCase
    [[nodiscard]] const LinearExtension& GetLeForConversion() const noexcept { return le_for_conversion_; }
    [[nodiscard]] const std::vector<std::uint64_t>& GetLabels() const noexcept { return labels_; }
    [[nodiscard]] const std::vector<std::uint64_t>& GetParent() const noexcept { return parent_; }
    [[nodiscard]] const std::vector<IdealSet>& GetIdeals() const noexcept { return ideals_; }
    // Il getter ora restituisce la referenza ai BitSet anziché al vector<vector>
    [[nodiscard]] const std::vector<BitSet>& GetImPred() const noexcept { return im_pred_; }
    [[nodiscard]] std::uint64_t GetRoot() const noexcept { return root_; }
    
    [[nodiscard]] std::vector<std::vector<std::uint64_t>> GetChildrenSortedLabel() const {
        std::vector<std::vector<std::uint64_t>> result(children_.size());
        
        for (size_t parent = 0; parent < children_.size(); ++parent) {
            const auto& childrenOf = children_[parent];
            auto& childrenOfResult = result[parent];
            childrenOfResult.reserve(childrenOf.size());
            
            std::vector<std::pair<std::uint64_t, std::uint64_t>> label_children;
            label_children.reserve(childrenOf.size());
            for (auto child : childrenOf) {
                label_children.emplace_back(labels_[child], child);
            }
            
            std::ranges::sort(label_children, std::greater<>{}, &std::pair<std::uint64_t, std::uint64_t>::first);
            
            for (const auto& [label, child] : label_children) {
                childrenOfResult.push_back(child);
            }
        }
        return result;
    }
    
private:
    std::uint64_t Left(std::int64_t n, IdealSet ideal) {
        std::uint64_t root = children_.size();
        children_.emplace_back();
        parent_.push_back(0);
        labels_.push_back(NO_LABEL);
        ideals_.push_back(ideal);
        
        if (n < 0) return root;
        
        IdealSet sub_ideal = ideal;
        sub_ideal.Unset(static_cast<std::uint64_t>(n));
        
        std::uint64_t r = Left(n - 1, std::move(sub_ideal));
        Right(static_cast<std::uint64_t>(n), r, root);
        
        parent_[r] = root;
        labels_[r] = static_cast<std::uint64_t>(n);
        children_[root].push_back(r);
        
        ideals_[r].Unset(static_cast<std::uint64_t>(n));
        
        return root;
    }
    
    void Right(std::uint64_t n, std::uint64_t r, std::uint64_t root) {
        for (auto child : children_[r]) {
            std::uint64_t label_s = labels_[child];
            
            // Check O(1) tramite BitSet nativo in im_pred_[n]
            if (!im_pred_[n].test(label_s)) {
                IdealSet sub_ideal = ideals_[root];
                sub_ideal.Unset(label_s);
                
                std::uint64_t t = children_.size();
                children_.emplace_back();
                parent_.push_back(0);
                labels_.push_back(NO_LABEL);
                ideals_.push_back(std::move(sub_ideal));
                
                parent_[t] = root;
                labels_[t] = label_s;
                children_[root].push_back(t);
                
                ideals_[t].Unset(label_s);
                Right(n, child, t);
            }
        }
    }
};
