#pragma once

#include "poset.h"
#include "tensor.h"
#include <vector>
#include <string>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <string_view>
#include <cstdint>

class BinaryVariablePOSet final : public POSet {
private:
    std::vector<std::string> variables_;
    std::uint64_t numero_profili_;
    std::vector<std::string> element_names_;
    
    // Costruttore privato (Pattern Factory)
    BinaryVariablePOSet() : POSet() {}
    
public:
    // -------------------------------------------------------------------------
    // Factory e Attributi Base
    // -------------------------------------------------------------------------
    [[nodiscard]] static std::unique_ptr<BinaryVariablePOSet> Build(const std::vector<std::string>& variables);
    [[nodiscard]] std::uint64_t NumberOfVariables() const noexcept;
    
    // -------------------------------------------------------------------------
    // Overrides Interfaccia POSet
    // -------------------------------------------------------------------------
    
    [[nodiscard]] std::uint64_t size() const noexcept override;
    
    [[nodiscard]] std::uint64_t GetElementId(std::string_view name) const override;
    [[nodiscard]] std::string_view GetElementName(std::uint64_t idx) const override;
    
    // -------------------------------------------------------------------------
    // Motore Relazionale On-The-Fly
    // -------------------------------------------------------------------------
    [[nodiscard]] bool IsLessOrEqual(std::uint64_t a, std::uint64_t b) const noexcept override;
    [[nodiscard]] bool GreaterThan(std::uint64_t e1, std::uint64_t e2) const noexcept override;
    [[nodiscard]] std::vector<std::pair<std::uint64_t, std::uint64_t>> Comparabilities() const override;
    [[nodiscard]] Tensor<std::uint8_t, 2> IncidenceMatrix() const override;
    [[nodiscard]] std::shared_ptr<POSet> Clone() const override;
    [[nodiscard]] std::unique_ptr<POSet> Dual() const override;
    [[nodiscard]] bool IsTotalOrder() const override;
    [[nodiscard]] Tensor<std::uint8_t, 2> CoverMatrix() const override;
    [[nodiscard]] bool IsExtensionOf(const POSet& p) const override;
    [[nodiscard]] BitSet DownSet(const std::vector<std::uint64_t>& els) const override;
    [[nodiscard]] bool IsDownSet(const std::vector<std::uint64_t>& els) const override;
    [[nodiscard]] BitSet UpSet(const std::vector<std::uint64_t>& els) const override;
    [[nodiscard]] const POSet::DATASTORE& UpSets() const override;
    [[nodiscard]] bool IsUpSet(const std::vector<std::uint64_t>& els) const override;
    [[nodiscard]] std::vector<std::pair<std::string, std::string>> OrderRelation() const override;
    [[nodiscard]] std::vector<std::pair<std::uint64_t, std::uint64_t>> CoverRelation() const override;
    [[nodiscard]] BitSet ComparabilitySetOf(std::uint64_t e) const override;
    [[nodiscard]] BitSet IncomparabilitySetOf(std::uint64_t e) const override;
    [[nodiscard]] BitSet Maximals() const override;
    [[nodiscard]] BitSet Minimals() const override;
    [[nodiscard]] bool IsMaximal(std::uint64_t e) const override;
    [[nodiscard]] bool IsMinimal(std::uint64_t e) const override;
    [[nodiscard]] std::vector<std::pair<std::uint64_t, std::uint64_t>> Incomparabilities() const override;
    void FirstLE(LinearExtension& le) const override;
    [[nodiscard]] bool IsComparable(std::uint64_t a, std::uint64_t b) const noexcept override;
    [[nodiscard]] Tensor<double, 2> BLSDominanceAbsolute() const override;
    [[nodiscard]] Tensor<double, 2> BLSDominanceRelative() const override;
    [[nodiscard]] std::unique_ptr<LinearExtensionGenerator> CreateLinearExtensionGenerator() override;
    // non implementati
    [[nodiscard]] LatticeOfIdeals* GetLatticeOfIdeals() override;
    [[nodiscard]] TreeOfIdeals* GetTreeOfIdeals() override;
    [[nodiscard]] std::string to_string(char delimiter = ';') const override;
      
};
