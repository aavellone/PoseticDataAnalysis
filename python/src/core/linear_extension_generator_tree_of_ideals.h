/**
 * @file linear_extension_generator_tree_of_ideals.h
 * @brief Header file for the Tree of Ideals linear extension generator.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "linear_extension_generator.h"

// Forward declaration to minimize compile-time dependencies.
class LatticeOfIdeals;

/**
 * @class LEGTreeOfIdeals
 * @brief Exact generator based on the Lattice of Ideals.
 *
 * @details Exploits the isomorphism between linear extensions and maximal paths
 * in the lattice of order ideals, derived from Birkhoff's representation theorem.
 */
class LEGTreeOfIdeals final : public LinearExtensionGenerator {
public:
    LEGTreeOfIdeals(std::uint64_t size, const LatticeOfIdeals& lattice);
    
    ~LEGTreeOfIdeals() override = default;
    
    void Start(std::uint64_t id = 0) override;
    
    void Next() override;
    
    [[nodiscard]] bool HasNext() override;
    
    [[nodiscard]] std::uint64_t NumberOfLe() const noexcept override {
        return 0;
    }
    
private:
    const LatticeOfIdeals& lattice_of_ideals_;      ///< Reference to the underlying lattice structure.
    
    std::vector<std::uint64_t> lattice_of_ideals_crossing_; ///< Tracks the path crossing within the lattice.
    
    std::vector<std::uint8_t> more_crossing_;
    
    std::vector<std::uint64_t> out_result_buffer_;
};
