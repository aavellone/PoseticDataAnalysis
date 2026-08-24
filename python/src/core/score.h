#pragma once

// MSVC: explicit standard includes (not pulled in transitively).
#include <utility>
#include <cstdint>
#include <cstddef>
#include <functional>


#include <vector>
#include <string>
#include <string_view>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <ranges>

class Score {
protected:
    // Hashing trasparente C++20 per permettere lookup con std::string_view
    // all'interno di una mappa che ha std::string_view come chiavi.
    struct StringHash {
        using is_transparent = void;
        [[nodiscard]] std::size_t operator()(std::string_view txt) const noexcept {
            return std::hash<std::string_view>{}(txt);
        }
    };
private:
    // Structure of Arrays (SoA): Ottimale per l'accesso contiguo in memoria (HPC)
    std::vector<std::string> names_;
    std::vector<double> scores_;
    
    // std::unordered_map è in media O(1) contro l'O(log N) di std::map.
    // Usiamo std::string_view per non allocare copie extra delle stringhe.
    std::unordered_map<std::string_view, std::uint64_t, StringHash, std::equal_to<>> name_to_eid_;
    
public:
    // Accetta genericamente std::map, std::unordered_map, o un vector di pairs
    template<typename Iterable>
    explicit Score(const Iterable& scores) {
        // 1. Pre-alloca per evitare riallocazioni
        std::vector<std::pair<std::string, double>> temp;
        temp.reserve(scores.size());
        
        for (const auto& [name, score] : scores) {
            temp.emplace_back(name, score);
        }
        
        // 2. Ordinamento O(N log N) idiomatico C++20 (decrescente per score)
        std::ranges::sort(temp, std::greater<>{}, &std::pair<std::string, double>::second);
        
        names_.reserve(temp.size());
        scores_.reserve(temp.size());
        name_to_eid_.reserve(temp.size());
        
        // 3. Popola i vettori SoA. Muoviamo le stringhe per evitare copie.
        for (auto& [name, score] : temp) {
            names_.push_back(std::move(name));
            scores_.push_back(score);
        }
        
        // 4. Popola la hash map.
        // ATTENZIONE: siccome names_ ha già la capacità corretta, non riallocherà.
        // Questo garantisce che gli std::string_view puntino a memoria stabile.
        for (std::uint64_t i = 0; i < names_.size(); ++i) {
            name_to_eid_.emplace(names_[i], i);
        }
    }
    
    // Usa [[nodiscard]] e noexcept per ottimizzazioni del compilatore
    [[nodiscard]] std::uint64_t size() const noexcept {
        return names_.size();
    }
    
    // Passa le stringhe by-view per zero-allocation
    [[nodiscard]] std::uint64_t getEID(std::string_view name) const {
        if (auto it = name_to_eid_.find(name); it != name_to_eid_.end()) {
            return it->second;
        }
        // Nota: usa std::runtime_error o la tua custom `MyException`
        throw std::runtime_error(std::string("Score error getEID: ") + std::string(name));
    }
    
    // Ritorna const reference invece che copia
    [[nodiscard]] const std::string& getName(std::uint64_t eid) const {
        if (eid >= names_.size()) { // BUG RISOLTO: prima era >= eid
            throw std::out_of_range("Score error getName: " + std::to_string(eid));
        }
        return names_[eid];
    }
    
    [[nodiscard]] double getScore(std::string_view name) const {
        // Usa getEID internamente, che è già validato
        return scores_[getEID(name)];
    }
    
    [[nodiscard]] double getScore(std::uint64_t eid) const {
        if (eid >= scores_.size()) {
            throw std::out_of_range("Score error getScore: " + std::to_string(eid));
        }
        return scores_[eid];
    }
};

