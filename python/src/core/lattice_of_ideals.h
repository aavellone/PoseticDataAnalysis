#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <cstdint>
#include <cstddef>

#include "my_exception.h"
#include "tree_of_ideals.h"
#include "bit_set.h"

#include <vector>
#include <queue>
#include <stack>
#include <optional>
#include <algorithm>
#include <ranges>
#include <stdexcept>


class LatticeOfIdeals {
private:
    std::uint64_t top_;
    std::uint64_t bottom_;
    std::vector<std::vector<std::uint64_t>> immediate_predecessors_;
    std::vector<std::vector<std::uint64_t>> immediate_successors_;
    
    std::vector<std::optional<std::uint64_t>> number_le_filter_;
    std::vector<std::optional<std::uint64_t>> number_le_ideal_;
    
    std::vector<std::vector<std::uint64_t>> with_label_;
    
    // HPC: Riferimento costante all'albero.
    // È safe perché POSet garantisce che l'albero sopravviva al reticolo.
    const TreeOfIdeals& tree_of_ideals_;
    
    // Metodo helper per inserimento unico rapido su vettori piccoli
    inline void InsertUnique(std::vector<std::uint64_t>& vec, std::uint64_t val) {
        if (std::ranges::find(vec, val) == vec.end()) {
            vec.push_back(val);
        }
    }
    
public:
    LatticeOfIdeals(const TreeOfIdeals& treeOfIdeals, std::uint64_t root)
    : tree_of_ideals_(treeOfIdeals) {
        
        auto childrenSortedLabel = tree_of_ideals_.GetChildrenSortedLabel();
        std::uint64_t num_nodes = tree_of_ideals_.GetLabels().size();
        
        top_ = root;
        immediate_predecessors_.assign(num_nodes, {});
        immediate_successors_.assign(num_nodes, {});
        
        for (auto pred : childrenSortedLabel[root]) {
            immediate_predecessors_[root].push_back(pred);
            InsertUnique(immediate_successors_[pred], root);
        }
        
        ComputeSort();
        
        for (std::int64_t p_with_label = static_cast<std::int64_t>(with_label_.size()) - 1; p_with_label >= 0; --p_with_label) {
            auto& with_label_k = with_label_[static_cast<std::size_t>(p_with_label)];
            for (auto i : with_label_k) {
                immediate_predecessors_[i].clear();
                std::uint64_t pi = tree_of_ideals_.GetParent()[i];
                
                std::uint64_t h = 0;
                std::uint64_t j_primo = immediate_predecessors_[pi][h];
                while (j_primo != i) {
                    std::uint64_t j = childrenSortedLabel[j_primo].front();
                    immediate_predecessors_[i].push_back(j);
                    InsertUnique(immediate_successors_[j], i);
                    j_primo = immediate_predecessors_[pi][++h];
                }
                for (auto pred : childrenSortedLabel[i]) {
                    immediate_predecessors_[i].push_back(pred);
                    InsertUnique(immediate_successors_[pred], i);
                }
                if (immediate_predecessors_[i].empty()) {
                    bottom_ = i;
                }
            }
            
            for (auto i : with_label_k) {
                std::uint64_t pi = tree_of_ideals_.GetParent()[i];
                // HPC: std::erase funziona perfettamente da C++20 sui vettori
                std::erase(childrenSortedLabel[pi], i);
            }
        }
    }
    
    void GetFromPath(const std::vector<std::uint64_t>& path,
                     std::vector<std::uint8_t>& more,
                     std::vector<std::uint64_t>& out_result) const {
        
        if (out_result.size() != path.size() || more.size() != path.size()) {
            throw MyException("Error: vector size mismatch.");
        }
        
        std::uint64_t val_start = 0;
        
        for (std::uint64_t k = 0; k < path.size(); ++k) {
            const auto& set_start = immediate_predecessors_[val_start];
            std::uint64_t val_end = set_start[path[k]];
            
            more[k] = (path[k] + 1 < set_start.size());
            
            const auto& ideal1 = tree_of_ideals_.GetIdeals()[val_start];
            const auto& ideal2 = tree_of_ideals_.GetIdeals()[val_end];
            
            std::uint64_t label = BitSet::FindFirstDifference(ideal1, ideal2);
            
            if (label == static_cast<std::uint64_t>(-1)) {
                throw MyException("Error: Corrupted path.");
            }
            
            // Richiama il getter della LinearExtension
            out_result[path.size() - k - 1] = tree_of_ideals_.GetLeForConversion().GetVal(label);
            
            val_start = val_end;
        }
    }
    
    void ComputeSort() {
        const auto& labels = tree_of_ideals_.GetLabels();
        with_label_.assign(tree_of_ideals_.GetImPred().size(), {});
        
        for (std::uint64_t lp = 0; lp < labels.size(); ++lp) {
            auto label = labels[lp];
            if (label != TreeOfIdeals::NO_LABEL) {
                InsertUnique(with_label_[label], lp);
            }
        }
    }
    
    std::uint64_t ComputeFiltersCount() { return AssignTopDown(); }
    void ComputeIdealsCount() { AssignBottomUp(); }
    
    [[nodiscard]] const auto& GetImmediatePredecessors() const noexcept { return immediate_predecessors_; }
    [[nodiscard]] const auto& GetImmediateSuccessors() const noexcept { return immediate_successors_; }
    [[nodiscard]] std::uint64_t GetBottom() const noexcept { return bottom_; }
    [[nodiscard]] std::uint64_t GetTop() const noexcept { return top_; }
    [[nodiscard]] const TreeOfIdeals& GetTreeOfIdeals() const noexcept { return tree_of_ideals_; }
    
    [[nodiscard]] const auto& GetFiltersCount() const noexcept { return number_le_filter_; }
    [[nodiscard]] const auto& GetIdealsCount() const noexcept { return number_le_ideal_; }
    
private:
    // HPC BUGFIX: Sostituita la BFS con l'Algoritmo di Kahn (Topological Sort).
    // Su un DAG, un nodo può essere sommato SOLO quando tutti i suoi padri
    // hanno completato i calcoli. La BFS causava doppie computazioni e calcoli parziali.
    void AssignBottomUp() {
        std::uint64_t num_nodes = immediate_predecessors_.size();
        number_le_ideal_.assign(num_nodes, std::nullopt);
        number_le_ideal_[bottom_] = 1;
        
        // 1. Calcoliamo il grado entrante (quanti padri ha ogni nodo in questa direzione)
        std::vector<std::uint64_t> in_degree(num_nodes, 0);
        for (std::uint64_t i = 0; i < num_nodes; ++i) {
            for (auto succ : immediate_successors_[i]) {
                in_degree[succ]++;
            }
        }
        
        std::queue<std::uint64_t> Q;
        Q.push(bottom_);
        
        // 2. Diffusione Topologica
        while (!Q.empty()) {
            auto v = Q.front();
            Q.pop();
            
            std::uint64_t val_v = number_le_ideal_[v].value_or(0);
            
            for (auto v_primo : immediate_successors_[v]) {
                std::uint64_t current_val = number_le_ideal_[v_primo].value_or(0);
                number_le_ideal_[v_primo] = current_val + val_v;
                
                // Diminuiamo il grado. Se arriva a 0, tutti i predecessori lo hanno aggiornato.
                if (--in_degree[v_primo] == 0) {
                    Q.push(v_primo);
                }
            }
        }
    }
    
    // AssignTopDown usava già un approccio DFS ricorsivo-iterativo, che per il conteggio top-down va bene.
    std::uint64_t AssignTopDown() {
        std::uint64_t num_nodes = immediate_predecessors_.size();
        number_le_filter_.assign(num_nodes, std::nullopt);
        
        number_le_filter_[top_] = 1;
        
        std::stack<std::uint64_t> S;
        S.push(bottom_);
        
        while (!S.empty()) {
            auto v = S.top();
            std::uint64_t e = 0;
            bool number_le_filter_for_all_imsuc = true;
            
            for (auto v_primo : immediate_successors_[v]) {
                if (number_le_filter_[v_primo].has_value()) {
                    e += number_le_filter_[v_primo].value();
                } else {
                    number_le_filter_for_all_imsuc = false;
                    S.push(v_primo);
                }
            }
            
            if (number_le_filter_for_all_imsuc) {
                if (v != top_) {
                    number_le_filter_[v] = e;
                }
                S.pop();
            }
        }
        
        return number_le_filter_[bottom_].value_or(0);
    }
};
