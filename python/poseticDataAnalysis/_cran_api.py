"""CRAN-compatible functional API.

The published Python package exposes exactly the functions of the CRAN R
package `poseticDataAnalysis` (https://cran.r-project.org/package=poseticDataAnalysis),
as module-level functions with the same names. Internally these delegate to the
POSet methods and to the generator/evaluation objects; the richer Python-only
surface (extra functions, method-style access, result classes) stays available
but is not part of the public CRAN-mirroring API.
"""

from __future__ import annotations

from typing import List, Optional

from ._evaluation import LEBubleyDyer

__all__ = [
    # structure / relation queries (free-function form, like R)
    "POSetElements", "DominanceMatrix", "CoverMatrix", "OrderRelation",
    "CoverRelation", "IncomparabilityRelation",
    "Dominates", "IsDominatedBy", "IsComparableWith", "IsIncomparableWith",
    "UpsetOf", "DownsetOf", "IsUpset", "IsDownset",
    "ComparabilitySetOf", "IncomparabilitySetOf",
    "POSetMaximals", "POSetMinimals", "IsMaximal", "IsMinimal",
    "POSetMeet", "POSetJoin", "IsExtensionOf",
    # generator / evaluation drivers
    "LEGet", "BubleyDyerMRP", "BubleyDyerEvaluation",
]


# --- structure -------------------------------------------------------------

def POSetElements(poset):
    """The element names of the poset."""
    return poset.elements


def DominanceMatrix(poset):
    """The dominance (incidence) matrix."""
    return poset.incidence_matrix()


def CoverMatrix(poset):
    """The cover (Hasse) matrix."""
    return poset.cover_matrix()


def OrderRelation(poset):
    """All order pairs (u, v) with u <= v."""
    return poset.order_relation()


def CoverRelation(poset):
    """Hasse-diagram covering pairs."""
    return poset.cover_relation()


def IncomparabilityRelation(poset):
    """All incomparable element pairs."""
    return poset.incomparabilities()


# --- pairwise comparisons --------------------------------------------------

def Dominates(poset, v1, v2):
    """Element-wise v1 >= v2."""
    return poset.dominates(v1, v2)


def IsDominatedBy(poset, v1, v2):
    """Element-wise v1 <= v2."""
    return poset.is_dominated_by(v1, v2)


def IsComparableWith(poset, v1, v2):
    """Element-wise comparability of v1 and v2."""
    return poset.is_comparable_with(v1, v2)


def IsIncomparableWith(poset, v1, v2):
    """Element-wise incomparability of v1 and v2."""
    return poset.is_incomparable_with(v1, v2)


# --- up / down sets --------------------------------------------------------

def UpsetOf(poset, elements):
    return poset.upset_of(elements)


def DownsetOf(poset, elements):
    return poset.downset_of(elements)


def IsUpset(poset, elements):
    return poset.is_upset(elements)


def IsDownset(poset, elements):
    return poset.is_downset(elements)


# --- per-element (in)comparability -----------------------------------------

def ComparabilitySetOf(poset, element):
    return poset.comparability_set_of(element)


def IncomparabilitySetOf(poset, element):
    return poset.incomparability_set_of(element)


# --- extremal --------------------------------------------------------------

def POSetMaximals(poset):
    return poset.maximals()


def POSetMinimals(poset):
    return poset.minimals()


def IsMaximal(poset, element):
    return poset.is_maximal(element)


def IsMinimal(poset, element):
    return poset.is_minimal(element)


def POSetMeet(poset, elements):
    """Greatest lower bound, or None."""
    return poset.meet(elements)


def POSetJoin(poset, elements):
    """Least upper bound, or None."""
    return poset.join(elements)


def IsExtensionOf(poset, other):
    return poset.is_extension_of(other)


# --- generator / evaluation drivers ----------------------------------------

def LEGet(generator, from_start: bool = True, n: Optional[int] = None,
          error: Optional[float] = None,
          output_every_sec: Optional[int] = None) -> List[List[str]]:
    """Draw linear extensions from a generator (exact or Bubley-Dyer).

    ``error`` applies to Bubley-Dyer generators only.
    """
    if isinstance(generator, LEBubleyDyer):
        return generator.get(from_start=from_start, n=n, error=error,
                             output_every=output_every_sec)
    return generator.get(from_start=from_start, n=n, output_every=output_every_sec)


def BubleyDyerMRP(generator, n: Optional[int] = None, error: Optional[float] = None,
                  output_every_sec: Optional[int] = None):
    """Run (or continue) a Bubley-Dyer MRP generator."""
    return generator.run(n=n, error=error, output_every=output_every_sec)


def BubleyDyerEvaluation(generator, n: Optional[int] = None, error: Optional[float] = None,
                         output_every_sec: Optional[int] = None):
    """Run (or continue) a Bubley-Dyer evaluation generator."""
    return generator.run(n=n, error=error, output_every=output_every_sec)
