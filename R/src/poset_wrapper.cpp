/**
 * @file poset_wrapper.cpp
 * @brief Implementazione 100% unique_ptr e semantica di spostamento.
 *
 * @author Alessandro Avellone
 * @version 7.0 (HPC Edition - 100% unique_ptr)
 */

#include "poset_wrapper.h"

#include <stdexcept>

#include "my_exception.h"
#include "poset.h"
#include "binary_variable_poset.h"
#include "lexicographic_product_poset.h"
#include "linear_poset.h"
#include "product_poset.h"
#include "bucket_poset.h"
#include "generic_functions.h"

const std::map<std::string, POSetWrap::FunctionLinearType> POSetWrap::kFunctionLinearMapType = {
    {"MutualRankingProbability", FunctionLinearType::kMutualRankingProbability},
    {"AverageHeight", FunctionLinearType::kAverageHeight},
    {"symmetric", FunctionLinearType::kSeparationSymmetric},
    {"asymmetricLower", FunctionLinearType::kSeparationAsymmetricLower},
    {"asymmetricUpper", FunctionLinearType::kSeparationAsymmetricUpper},
    {"RFunction", FunctionLinearType::kRFunction},
};

// Implementazione del distruttore nel .cpp
POSetWrap::~POSetWrap() = default;

// =============================================================================
// Factories Implementations (std::move enabled)
// =============================================================================

std::unique_ptr<POSetWrap> POSetWrap::BuildPoset(
                                                   const std::vector<std::string>& elements,
                                                   const std::vector<std::pair<std::string, std::string>>& comparabilities) {
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kPoset;
    wrapper->poset_ = POSet::Build(elements, comparabilities);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildLinear(const std::vector<std::string>& elements) {
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kLin;
    wrapper->poset_ = LinearPOSet::Build(elements);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildProduct(const std::vector<const POSetWrap*>& posets) {
    if (posets.size() < 2) throw MyException("Servono almeno due POSet per il Prodotto.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kProd;
    
    std::unique_ptr<POSet> current = ProductPOSet::Build(*posets[0]->GetPOSet(), *posets[1]->GetPOSet());
    for (std::size_t i = 2; i < posets.size(); ++i) {
        current = ProductPOSet::Build(*current, *posets[i]->GetPOSet());
    }
    
    wrapper->poset_ = std::move(current); // Move semantics
    return wrapper;
}


std::unique_ptr<POSetWrap> POSetWrap::BuildLexicographicProduct(const std::vector<const POSetWrap*>& posets) {
    if (posets.size() < 2) throw MyException("Servono almeno due POSet per LexicographicProduct.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kLexProd;
    
    std::unique_ptr<POSet> current = LexicographicProductPOSet::Build(*posets[0]->GetPOSet(), *posets[1]->GetPOSet());
    for (std::size_t i = 2; i < posets.size(); ++i) {
        current = LexicographicProductPOSet::Build(*current, *posets[i]->GetPOSet());
    }
    
    wrapper->poset_ = std::move(current);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildBinaryVariablePOSet(const std::vector<std::string>& elements) {
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kBinaryVariable;
    wrapper->poset_ = BinaryVariablePOSet::Build(elements);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildBucketPOSet(
                                                       const std::vector<std::string>& elements,
                                                       const std::vector<std::pair<std::string, std::string>>& buckets) {
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kBucket;
    wrapper->poset_ = BucketPOSet::Build(elements, buckets);
    return wrapper;
}


std::unique_ptr<POSetWrap> POSetWrap::BuildIntersection(const std::vector<const POSetWrap*>& posets) {
    if (posets.size() < 2) throw MyException("Servono almeno due POSet per l'Intersezione.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kIntersection;
    
    std::unique_ptr<POSet> current = POSet::Intersection(*posets[0]->GetPOSet(), *posets[1]->GetPOSet());
    for (std::size_t i = 2; i < posets.size(); ++i) {
        current = POSet::Intersection(*current, *posets[i]->GetPOSet());
    }
    
    wrapper->poset_ = std::move(current);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildLinearSum(const std::vector<const POSetWrap*>& posets) {
    if (posets.size() < 2) throw MyException("Servono almeno due POSet per la LinearSum.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kLinearSum;
    
    std::unique_ptr<POSet> current = POSet::LinearSum(*posets[0]->GetPOSet(), *posets[1]->GetPOSet());
    for (std::size_t i = 2; i < posets.size(); ++i) {
        current = POSet::LinearSum(*current, *posets[i]->GetPOSet());
    }
    
    wrapper->poset_ = std::move(current);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildDisjointSum(const std::vector<const POSetWrap*>& posets) {
    if (posets.size() < 2) throw MyException("Servono almeno due POSet per la DisjointSum.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kDisjointSum;
    
    std::unique_ptr<POSet> current = POSet::DisjointSum(*posets[0]->GetPOSet(), *posets[1]->GetPOSet());
    for (std::size_t i = 2; i < posets.size(); ++i) {
        current = POSet::DisjointSum(*current, *posets[i]->GetPOSet());
    }
    
    wrapper->poset_ = std::move(current);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildLiftingPOSet(const POSetWrap* poset, const std::string& new_element) {
    if (!poset || !poset->GetPOSet()) throw MyException("POSet non valido per il Lifting.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kLifting;
    wrapper->poset_ = poset->GetPOSet()->Lifting(new_element);
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildDualPOSet(const POSetWrap* poset) {
    if (!poset || !poset->GetPOSet()) throw MyException("POSet di base non valido per il Dual.");
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kDual;
    wrapper->poset_ = poset->GetPOSet()->Dual();
    return wrapper;
}

std::unique_ptr<POSetWrap> POSetWrap::BuildFencePOSet(const std::vector<std::string>& elements, bool orientation) {
    // Generazione delle comparabilità a zig-zag
    std::vector<std::pair<std::string, std::string>> comparabilities;
    
    if (elements.size() > 1) {
        comparabilities.reserve(elements.size() - 1);
        bool current_orientation = orientation;
        
        for (std::size_t i = 0; i < elements.size() - 1; ++i) {
            if (current_orientation) {
                // elements[i] < elements[i+1]
                comparabilities.emplace_back(elements[i], elements[i + 1]);
            } else {
                // elements[i] > elements[i+1]
                comparabilities.emplace_back(elements[i + 1], elements[i]);
            }
            // Alterna l'orientamento al prossimo passo
            current_orientation = !current_orientation;
        }
    }
    
    auto wrapper = std::unique_ptr<POSetWrap>(new POSetWrap());
    wrapper->type_ = PosetType::kFence;
    // Si affida al costruttore del POSet Generico passando la lista archi creata ad-hoc
    wrapper->poset_ = POSet::Build(elements, comparabilities);
    
    return wrapper;
}



