"""Separation, fuzzy measures, dimensionality reduction and FOD analysis.

Mirrors the corresponding R functions. Separation on a poset reuses the same
evaluation machinery (the R package does the same); the lexicographic, fuzzy,
dimensionality-reduction and First Order Dominance functions call dedicated
core routines.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Union

from . import _core
from ._poset import POSet
from ._evaluation import EvaluationResult, Seed, _seed_arg

__all__ = [
    "SEPARATION_METRICS",
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
    "OptimalBidimensionalEmbedding",
    "BidimentionalPosetRepresentation",
    "FirstOrderDominanceAnalysis",
]

SEPARATION_METRICS = ("symmetric", "asymmetricLower", "asymmetricUpper")
_FUZZY_CORE = ("symmetric", "asymmetricLower", "asymmetricUpper")
_FOD_METRIC_NAMES = ("Dominance", "MannWhitneyDominance",
                     "MannWhitneyInferentialDominance", "MRP")


def _check_sep(quali):
    q = list(quali)
    for name in q:
        if name not in SEPARATION_METRICS:
            raise ValueError(
                f"unknown separation type {name!r}; valid: {', '.join(SEPARATION_METRICS)}")
    return q


# ---------------------------------------------------------------------------
# Separation on a poset (reuses the evaluation engine)
# ---------------------------------------------------------------------------

def ExactSeparation(poset: POSet, quali: Sequence[str] = ("symmetric",),
                    output_every: Optional[int] = None) -> EvaluationResult:
    """Exact separation matrices over every linear extension."""
    return EvaluationResult(_core.exact_evaluation(poset.handle, _check_sep(quali), output_every))


class BuildBubleyDyerSeparationGenerator:
    """Stateful MCMC estimator of separation matrices (Bubley–Dyer)."""

    __slots__ = ("_handle", "_poset", "seed", "_quali")

    def __init__(self, poset: POSet, quali: Sequence[str] = ("symmetric",),
                 seed: Seed = None) -> None:
        self._poset = poset
        self._quali = _check_sep(quali)
        self._handle, self.seed = _core.build_bubley_dyer_evaluation_generator(
            poset.handle, self._quali, _seed_arg(seed))

    def run(self, n: Optional[int] = None, error: Optional[float] = None,
            output_every: Optional[int] = None) -> EvaluationResult:
        return EvaluationResult(
            _core.bubley_dyer_evaluation(self._handle, n, error, output_every))


def BubleyDyerSeparation(generator: BuildBubleyDyerSeparationGenerator,
                         n: Optional[int] = None, error: Optional[float] = None,
                         output_every: Optional[int] = None) -> EvaluationResult:
    """Run (or continue) a Bubley–Dyer separation generator."""
    return generator.run(n=n, error=error, output_every=output_every)


# ---------------------------------------------------------------------------
# Lexicographic separation / MRP
# ---------------------------------------------------------------------------

def _build_modalita(nvar: int, deg) -> List[List[str]]:
    nvar = int(nvar)
    if nvar <= 0:
        raise ValueError("'nvar' must be a positive integer.")
    if isinstance(deg, int):
        return [[str(j) for j in range(1, deg + 1)] for _ in range(nvar)]
    deg = list(deg)
    # a flat list of level names shared by all variables
    if deg and all(isinstance(x, str) for x in deg):
        return [list(deg) for _ in range(nvar)]
    # a per-variable list (of ints -> level counts, or of level-name lists)
    if len(deg) != nvar:
        raise ValueError(f"'deg' must have length {nvar}.")
    out = []
    for x in deg:
        if isinstance(x, int):
            out.append([str(j) for j in range(1, x + 1)])
        else:
            out.append([str(v) for v in x])
    return out


def LexSeparation(nvar: int, deg,
                  types: Sequence[str] = ("symmetric",)) -> Dict[str, object]:
    """Lexicographic separation for a product of ``nvar`` discrete variables.

    ``deg`` may be an int (levels per variable), a list of level names shared by
    all variables, or a per-variable list. Returns the requested ``types``
    (``symmetric``, ``asymmetricLower``, ``asymmetricUpper``, ``vertical``,
    ``horizontal``) plus ``labels``.
    """
    valid = ("symmetric", "asymmetricLower", "asymmetricUpper", "vertical", "horizontal")
    types = list(types)
    for t in types:
        if t not in valid:
            raise ValueError(f"unknown type {t!r}; valid: {', '.join(valid)}")
    modalita = _build_modalita(nvar, deg)
    full = _core.lex_separation(modalita)
    out = {t: full[t] for t in types}
    out["labels"] = full["labels"]
    return out


def LexMRP(nvar: int, deg) -> Dict[str, object]:
    """Mutual Ranking Probability matrix for a lexicographic product."""
    return _core.lex_mrp(_build_modalita(nvar, deg))


# ---------------------------------------------------------------------------
# Fuzzy separation / in-betweenness (on a dominance matrix)
# ---------------------------------------------------------------------------

def _check_fuzzy(quali):
    q = list(quali)
    for name in q:
        if name not in _FUZZY_CORE:
            raise ValueError(
                f"fuzzy type {name!r} not supported; valid: {', '.join(_FUZZY_CORE)} "
                "(derived 'vertical'/'horizontal' are not available in the Python binding).")
    return q


def FuzzySeparation(dom: Sequence[Sequence[float]], elements: Sequence[str],
                    norm: str = "minimum",
                    quali: Sequence[str] = ("symmetric",)) -> Dict[str, object]:
    """Fuzzy separation matrices for a dominance matrix (norm: 'minimum'/'product')."""
    return _core.fuzzy_separation([list(r) for r in dom], list(elements), norm, _check_fuzzy(quali))


def FuzzyInBetweenness(dom: Sequence[Sequence[float]], elements: Sequence[str],
                       norm: str = "minimum",
                       quali: Sequence[str] = ("symmetric",)) -> Dict[str, object]:
    """Fuzzy in-betweenness 3D arrays for a dominance matrix."""
    return _core.fuzzy_inbetweenness([list(r) for r in dom], list(elements), norm, _check_fuzzy(quali))


def FuzzySeparationMinMax(dom, elements, quali=("symmetric",)):
    """Fuzzy separation with the Gödel (min/max) norm."""
    return FuzzySeparation(dom, elements, "minimum", quali)


def FuzzySeparationProbabilistic(dom, elements, quali=("symmetric",)):
    """Fuzzy separation with the product / probabilistic-sum norm."""
    return FuzzySeparation(dom, elements, "product", quali)


def FuzzyInBetweennessMinMax(dom, elements, quali=("symmetric",)):
    """Fuzzy in-betweenness with the Gödel (min/max) norm."""
    return FuzzyInBetweenness(dom, elements, "minimum", quali)


def FuzzyInBetweennessProbabilistic(dom, elements, quali=("symmetric",)):
    """Fuzzy in-betweenness with the product / probabilistic-sum norm."""
    return FuzzyInBetweenness(dom, elements, "product", quali)


# ---------------------------------------------------------------------------
# Dimensionality reduction
# ---------------------------------------------------------------------------

def _lpom_int(lpom_strategy: str) -> int:
    if lpom_strategy == "absolute":
        return 0
    if lpom_strategy == "relative":
        return 1
    raise ValueError("lpom_strategy must be 'absolute' or 'relative'.")


def _drop_zero_weights(profile, weights):
    profile = [list(r) for r in profile]
    weights = list(weights)
    keep = [i for i, w in enumerate(weights) if w != 0]
    if len(keep) != len(weights):
        profile = [profile[i] for i in keep]
        weights = [weights[i] for i in keep]
    return profile, weights


def OptimalBidimensionalEmbedding(profile: Sequence[Sequence[int]], weights: Sequence[float],
                                  lpom_strategy: str = "absolute",
                                  output_every: Optional[int] = None,
                                  thread_share: float = 1.0,
                                  loss: str = "LB") -> Dict[str, object]:
    """Search all variable permutations for the best 2D embedding of the profiles."""
    if len(weights) != len(profile):
        raise ValueError("'weights' length must equal the number of profile rows.")
    profile, weights = _drop_zero_weights(profile, weights)
    return _core.dimensionality_reduction(profile, weights, loss, _lpom_int(lpom_strategy),
                                          output_every, float(thread_share))


def BidimentionalPosetRepresentation(profile: Sequence[Sequence[int]], weights: Sequence[float],
                                     variables_priority: Sequence[int],
                                     lpom_strategy: str = "absolute",
                                     loss: str = "LB") -> Dict[str, object]:
    """2D embedding of the profiles for a single, fixed variable priority (1-indexed)."""
    if len(weights) != len(profile):
        raise ValueError("'weights' length must equal the number of profile rows.")
    profile, weights = _drop_zero_weights(profile, weights)
    return _core.bidimensional_poset_representation(profile, weights, loss,
                                                    _lpom_int(lpom_strategy),
                                                    list(variables_priority))


# ---------------------------------------------------------------------------
# First Order Dominance analysis
# ---------------------------------------------------------------------------

def FirstOrderDominanceAnalysis(
        posets: Union[POSet, Sequence[POSet]],
        freq_matrix: Sequence[Sequence[float]],
        row_labels: Sequence[str],
        col_labels: Sequence[str],
        metrics: Sequence[str] = ("Dominance",),
        subpopulation_count: Optional[Sequence[float]] = None,
        total_bins: int = 0,
        count: Optional[int] = None,
        seed: Seed = None,
        output_every: Optional[int] = None,
        sep: str = "_",
        linear_extensions: Optional[Sequence[Sequence[str]]] = None,
        n_threads: Optional[int] = None) -> Dict[str, object]:
    """First Order Dominance analysis over one or more posets.

    ``posets`` is a single :class:`POSet` or a sequence of chains (for the
    lexicographic case). ``freq_matrix`` is a numeric matrix whose rows are
    ``row_labels`` (observed profiles, joined by ``sep`` in the multi-poset case)
    and whose columns are ``col_labels`` (groups). ``metrics`` selects among
    ``Dominance``, ``MannWhitneyDominance``, ``MannWhitneyInferentialDominance``,
    ``MRP``.

    Returns a dict with ``metrics`` (list of per-metric ``{name, fodClosed,
    fodMatrix, binMatrix, rows, cols}``), ``LEType``, ``le_count``, ``labels``
    and ``seed``.
    """
    handles = [posets.handle] if isinstance(posets, POSet) else [p.handle for p in posets]
    for name in metrics:
        if name not in _FOD_METRIC_NAMES:
            raise ValueError(
                f"unknown metric {name!r}; valid: {', '.join(_FOD_METRIC_NAMES)}")
    flags = [1 if name in metrics else 0 for name in _FOD_METRIC_NAMES]

    freq = [list(r) for r in freq_matrix]
    subpop = list(subpopulation_count) if subpopulation_count is not None else None
    les = ([list(le) for le in linear_extensions]
           if linear_extensions is not None else None)
    bins = int(total_bins) if total_bins else None

    return _core.first_order_dominance_analysis(
        handles, freq, list(row_labels), list(col_labels), flags, subpop, bins,
        count, _seed_arg(seed), output_every, sep, les, n_threads)
