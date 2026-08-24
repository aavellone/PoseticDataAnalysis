#include "rwrapper.h"

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>
#include <R.h>

extern "C" {
    static const R_CallMethodDef CallEntries[] = {
        {"BruggemannLercheSorensenDominance", (DL_FUNC) &BruggemannLercheSorensenDominance, 2},
        {"BubleyDyerEvaluation", (DL_FUNC) &BubleyDyerEvaluation, 4},
        {"BubleyDyerMRP", (DL_FUNC) &BubleyDyerMRP, 4},
        {"BubleyDyerSeparation", (DL_FUNC) &BubleyDyerSeparation, 4},
        {"BuildBinaryVariablePOSet", (DL_FUNC) &BuildBinaryVariablePOSet, 1},
        {"BuildBubleyDyerEvaluationGenerator", (DL_FUNC) &BuildBubleyDyerEvaluationGenerator, 3},
        {"BuildBubleyDyerLEGenerator", (DL_FUNC) &BuildBubleyDyerLEGenerator, 2},
        {"BuildBubleyDyerMRPGenerator", (DL_FUNC) &BuildBubleyDyerMRPGenerator, 2},
        {"BuildBubleyDyerSeparationGenerator", (DL_FUNC) &BuildBubleyDyerSeparationGenerator, 3},
        {"BuildBucketPOSet", (DL_FUNC) &BuildBucketPOSet, 2},
        {"BuildCrownPOSet", (DL_FUNC) &BuildCrownPOSet, 2},
        {"BuildDisjointSumPOSet", (DL_FUNC) &BuildDisjointSumPOSet, 1},
        {"BuildDualPOSet", (DL_FUNC) &BuildDualPOSet, 1},
        {"BuildFencePOSet", (DL_FUNC) &BuildFencePOSet, 2},
        {"BuildIntersectionPOSet", (DL_FUNC) &BuildIntersectionPOSet, 1},
        {"BuildLEGenerator", (DL_FUNC) &BuildLEGenerator, 1},
        {"BuildLexicographicProductPOSet", (DL_FUNC) &BuildLexicographicProductPOSet, 1},
        {"BuildLiftingPOSet", (DL_FUNC) &BuildLiftingPOSet, 2},
        {"BuildLinearPOSet", (DL_FUNC) &BuildLinearPOSet, 1},
        {"BuildLinearSumPOSet", (DL_FUNC) &BuildLinearSumPOSet, 1},
        {"BuildPOSet", (DL_FUNC) &BuildPOSet, 2},
        {"BuildProductPOSet", (DL_FUNC) &BuildProductPOSet, 1},
        {"ComparabilitySetOf", (DL_FUNC) &ComparabilitySetOf, 2},
        {"CoverMatrix", (DL_FUNC) &CoverMatrix, 1},
        {"CoverRelation", (DL_FUNC) &CoverRelation, 1},
        {"Dominates", (DL_FUNC) &Dominates, 3},
        {"DownsetOf", (DL_FUNC) &DownsetOf, 2},
        {"Elements", (DL_FUNC) &Elements, 1},
        {"ExactEvaluation", (DL_FUNC) &ExactEvaluation, 3},
        {"ExactMRP", (DL_FUNC) &ExactMRP, 2},
        {"ExactSeparation", (DL_FUNC) &ExactSeparation, 3},
        {"FirstOrderDominanceAnalysis", (DL_FUNC) &FirstOrderDominanceAnalysis, 11},
        {"FuzzyInBetweenness", (DL_FUNC) &FuzzyInBetweenness, 4},
        {"FuzzySeparation", (DL_FUNC) &FuzzySeparation, 4},
        {"IncidenceMatrix", (DL_FUNC) &IncidenceMatrix, 1},
        {"Incomparabilities", (DL_FUNC) &Incomparabilities, 1},
        {"IncomparabilitySetOf", (DL_FUNC) &IncomparabilitySetOf, 2},
        {"isAntisymmetric", (DL_FUNC) &isAntisymmetric, 1},
        {"IsComparableWith", (DL_FUNC) &IsComparableWith, 3},
        {"IsDominatedBy", (DL_FUNC) &IsDominatedBy, 3},
        {"IsDownset", (DL_FUNC) &IsDownset, 2},
        {"IsExtensionOf", (DL_FUNC) &IsExtensionOf, 2},
        {"IsIncomparableWith", (DL_FUNC) &IsIncomparableWith, 3},
        {"IsMaximal", (DL_FUNC) &IsMaximal, 2},
        {"IsMinimal", (DL_FUNC) &IsMinimal, 2},
        {"isPartialOrder", (DL_FUNC) &isPartialOrder, 2},
        {"isPreorder", (DL_FUNC) &isPreorder, 2},
        {"isReflexive", (DL_FUNC) &isReflexive, 2},
        {"isSymmetric", (DL_FUNC) &isSymmetric, 1},
        {"isTransitive", (DL_FUNC) &isTransitive, 1},
        {"IsUpset", (DL_FUNC) &IsUpset, 2},
        {"Join", (DL_FUNC) &Join, 2},
        {"LEGBubleyDyerGet", (DL_FUNC) &LEGBubleyDyerGet, 5},
        {"LEGGet", (DL_FUNC) &LEGGet, 4},
        {"Maximal", (DL_FUNC) &Maximal, 1},
        {"Meet", (DL_FUNC) &Meet, 2},
        {"Minimal", (DL_FUNC) &Minimal, 1},
        {"OrderRelation", (DL_FUNC) &OrderRelation, 1},
        {"ReflexiveClosure", (DL_FUNC) &ReflexiveClosure, 2},
        {"RunBidimentionalPosetRepresentation", (DL_FUNC) &RunBidimentionalPosetRepresentation, 5},
        {"RunDimensionalityReduction", (DL_FUNC) &RunDimensionalityReduction, 6},
        {"RunLexMRP", (DL_FUNC) &RunLexMRP, 1},
        {"RunLexSeparation", (DL_FUNC) &RunLexSeparation, 1},
        {"TransitiveClosure", (DL_FUNC) &TransitiveClosure, 1},
        {"UpsetOf", (DL_FUNC) &UpsetOf, 2},
        {NULL, NULL, 0}
    };

    void R_init_poseticDataAnalysis(DllInfo *dll) {
        R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);
        R_forceSymbols(dll, TRUE);
    }
}

