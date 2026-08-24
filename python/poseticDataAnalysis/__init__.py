"""poseticDataAnalysis — build and analyse partially ordered sets.

Native Python binding (CPython C-API, no third-party runtime dependencies) over
the same C++20 core used by the R package of the same name.

The public API mirrors the CRAN R package
(https://cran.r-project.org/package=poseticDataAnalysis): the module exposes the
same functions, with the same names. Additional Python-only capabilities exist
in the implementation but are intentionally not part of this public surface.
"""

from __future__ import annotations

# --- POSet constructors ----------------------------------------------------
from ._poset import (
    POSet,
    LinearPOSet,
    ProductPOSet,
    LexicographicProductPOSet,
    IntersectionPOSet,
    LinearSumPOSet,
    DisjointSumPOSet,
    LiftingPOSet,
    BinaryVariablePOSet,
    FencePOSet,
    CrownPOSet,
    DualPOSet,
)

# --- Relation property checks / closures -----------------------------------
from ._relations import (
    IsReflexive,
    IsSymmetric,
    IsAntisymmetric,
    IsTransitive,
    IsPreorder,
    IsPartialOrder,
    TransitiveClosure,
    ReflexiveClosure,
)

# --- Linear extensions / MRP / evaluation ----------------------------------
from ._evaluation import (
    LEGenerator,
    LEBubleyDyer,
    ExactMRP,
    BubleyDyerMRPGenerator,
    ExactEvaluation,
    BuildBubleyDyerEvaluationGenerator,
    BLSDominance,
)

# --- Separation / fuzzy / lexicographic / dimensionality reduction ---------
from ._advanced import (
    ExactSeparation,
    BuildBubleyDyerSeparationGenerator,
    BubleyDyerSeparation,
    LexSeparation,
    LexMRP,
    FuzzySeparation,
    FuzzyInBetweenness,
    FuzzySeparationMinMax,
    FuzzySeparationProbabilistic,
    FuzzyInBetweennessMinMax,
    FuzzyInBetweennessProbabilistic,
    OptimalBidimensionalEmbedding,
    BidimentionalPosetRepresentation,
)

# --- Free-function query / relation / generator API (CRAN parity) ----------
from ._cran_api import (
    POSetElements,
    DominanceMatrix,
    CoverMatrix,
    OrderRelation,
    CoverRelation,
    IncomparabilityRelation,
    Dominates,
    IsDominatedBy,
    IsComparableWith,
    IsIncomparableWith,
    UpsetOf,
    DownsetOf,
    IsUpset,
    IsDownset,
    ComparabilitySetOf,
    IncomparabilitySetOf,
    POSetMaximals,
    POSetMinimals,
    IsMaximal,
    IsMinimal,
    POSetMeet,
    POSetJoin,
    IsExtensionOf,
    LEGet,
    BubleyDyerMRP,
    BubleyDyerEvaluation,
)

__version__ = "1.1.1"

# Exactly the exported functions of the CRAN R package (plus __version__).
__all__ = [
    # POSet constructors
    "POSet",
    "LinearPOSet",
    "ProductPOSet",
    "LexicographicProductPOSet",
    "IntersectionPOSet",
    "LinearSumPOSet",
    "DisjointSumPOSet",
    "LiftingPOSet",
    "BinaryVariablePOSet",
    "FencePOSet",
    "CrownPOSet",
    "DualPOSet",
    # relation checks / closures
    "IsReflexive",
    "IsSymmetric",
    "IsAntisymmetric",
    "IsTransitive",
    "IsPreorder",
    "IsPartialOrder",
    "TransitiveClosure",
    "ReflexiveClosure",
    # structure / relation queries
    "POSetElements",
    "DominanceMatrix",
    "CoverMatrix",
    "OrderRelation",
    "CoverRelation",
    "IncomparabilityRelation",
    "Dominates",
    "IsDominatedBy",
    "IsComparableWith",
    "IsIncomparableWith",
    "UpsetOf",
    "DownsetOf",
    "IsUpset",
    "IsDownset",
    "ComparabilitySetOf",
    "IncomparabilitySetOf",
    "POSetMaximals",
    "POSetMinimals",
    "IsMaximal",
    "IsMinimal",
    "POSetMeet",
    "POSetJoin",
    "IsExtensionOf",
    # linear extensions / MRP / evaluation
    "LEGenerator",
    "LEBubleyDyer",
    "LEGet",
    "ExactMRP",
    "BubleyDyerMRPGenerator",
    "BubleyDyerMRP",
    "ExactEvaluation",
    "BuildBubleyDyerEvaluationGenerator",
    "BubleyDyerEvaluation",
    "BLSDominance",
    # separation / fuzzy / lexicographic
    "ExactSeparation",
    "BuildBubleyDyerSeparationGenerator",
    "BubleyDyerSeparation",
    "LexSeparation",
    "LexMRP",
    "FuzzySeparation",
    "FuzzyInBetweenness",
    "FuzzySeparationMinMax",
    "FuzzySeparationProbabilistic",
    "FuzzyInBetweennessMinMax",
    "FuzzyInBetweennessProbabilistic",
    # dimensionality reduction
    "OptimalBidimensionalEmbedding",
    "BidimentionalPosetRepresentation",
    # metadata
    "__version__",
]
