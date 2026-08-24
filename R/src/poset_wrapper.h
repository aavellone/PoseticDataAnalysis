/**
 * @file poset_wrapper.h
 * @brief Classe wrapper HPC per la gestione ottimizzata dei POSet.
 *
 * @details Fornisce un'interfaccia basata sul pattern Factory Method.
 * Gestisce l'ownership di POSet ESCLUSIVAMENTE tramite std::unique_ptr.
 * Zero-Cost Abstractions: elimina ogni overhead di atomic reference counting.
 *
 * @author Alessandro Avellone
 * @version 7.0 (HPC Edition - 100% unique_ptr)
 */

#pragma once

#include "poset.h"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/**
 * @class POSetWrap
 * @brief Wrapper sicuro e performante delle istanze POSet.
 */
class POSetWrap {
  public:
    enum class PosetType {
        kPoset,
        kLin,
        kProd,
        kLexProd,
        kBucket,
        kBinaryVariable,
        kIntersection,
        kLinearSum,
        kDisjointSum,
        kLifting,
        kFence,
        kDual
    };

    enum class FunctionLinearType {
        kMutualRankingProbability,
        kAverageHeight,
        kSeparationAsymmetricLower,
        kSeparationAsymmetricUpper,
        kSeparationSymmetric,
        kDominance,
        kMannWhitneyDominance,
        kMannWhitneyInferentialDominance,
        kRFunction,
    };

    static const std::map<std::string, FunctionLinearType>
        kFunctionLinearMapType;

    ~POSetWrap();

    // =====================================================================
    // Factories - Solo unique_ptr e const reference
    // =====================================================================

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildPoset(
        const std::vector<std::string> &elements,
        const std::vector<std::pair<std::string, std::string>>
            &comparabilities);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildBucketPOSet(
        const std::vector<std::string> &elements,
        const std::vector<std::pair<std::string, std::string>> &buckets);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildLinear(
        const std::vector<std::string> &elements);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildProduct(
        const std::vector<const POSetWrap *> &posets);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildLexicographicProduct(
        const std::vector<const POSetWrap *> &posets);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildBinaryVariablePOSet(
        const std::vector<std::string> &elements);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildIntersection(
        const std::vector<const POSetWrap *> &posets);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildLinearSum(
        const std::vector<const POSetWrap *> &posets);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildDisjointSum(
        const std::vector<const POSetWrap *> &posets);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildLiftingPOSet(
        const POSetWrap *poset, const std::string &new_element);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildFencePOSet(
        const std::vector<std::string> &elements, bool orientation);

    [[nodiscard]] static std::unique_ptr<POSetWrap> BuildDualPOSet(
        const POSetWrap *poset);

    // =====================================================================
    // Accessors
    // =====================================================================

    [[nodiscard]] POSet *GetPOSet() const noexcept { return poset_.get(); }
    [[nodiscard]] std::uint64_t GetId() const noexcept { return id_; }
    [[nodiscard]] PosetType GetType() const noexcept { return type_; }

  private:
    POSetWrap() = default;

    std::uint64_t id_ = 0;
    PosetType type_ = PosetType::kPoset;

    // Ownership pura: il wrapper possiede il POSet univocamente.
    std::unique_ptr<POSet> poset_;
};
