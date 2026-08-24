"""Linear-extension generators and evaluation over posets.

Mirrors the R package's generator/evaluation API. Long computations use the
same C++20 core; 64-bit seeds are carried as decimal strings (Python ints of
any size are accepted and converted losslessly).

Each generator object keeps a reference to its source :class:`POSet` so the
native POSet is never freed while the generator still points at it (the analog
of ``R_SetExternalPtrProtected`` in the R package).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Union

from . import _core
from ._poset import POSet

__all__ = [
    "LEGenerator",
    "LEBubleyDyer",
    "ExactMRP",
    "BubleyDyerMRPGenerator",
    "ExactEvaluation",
    "BuildBubleyDyerEvaluationGenerator",
    "BLSDominance",
    "MRPResult",
    "LabeledMatrix",
    "EvaluationResult",
    "VALID_METRICS",
]

VALID_METRICS = (
    "MutualRankingProbability",
    "AverageHeight",
    "symmetric",
    "asymmetricLower",
    "asymmetricUpper",
)

Seed = Optional[Union[int, str]]


def _seed_arg(seed: Seed) -> object:
    """Normalise a seed for the native layer (int stays int, str stays str)."""
    return seed  # None | int | str are all accepted by the C layer


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

class MRPResult:
    """Mutual Ranking Probability result: a square matrix over the elements."""

    __slots__ = ("matrix", "elements", "n")

    def __init__(self, d: dict) -> None:
        self.matrix: List[List[float]] = d["mrp"]
        self.elements: List[str] = d["elements"]
        self.n: int = d["n"]

    def __repr__(self) -> str:
        return f"MRPResult(n_extensions={self.n}, size={len(self.elements)})"


class LabeledMatrix:
    """A result matrix with row and column labels."""

    __slots__ = ("name", "matrix", "rows", "cols")

    def __init__(self, d: dict) -> None:
        self.name: str = d["name"]
        self.matrix: List[List[float]] = d["matrix"]
        self.rows: List[str] = d["rows"]
        self.cols: List[str] = d["cols"]

    def __repr__(self) -> str:
        return f"LabeledMatrix(name={self.name!r}, shape=({len(self.rows)}, {len(self.cols)}))"


class EvaluationResult:
    """Result of a multi-metric evaluation: metrics by name, plus n."""

    __slots__ = ("metrics", "n", "_order")

    def __init__(self, d: dict) -> None:
        self.metrics: Dict[str, LabeledMatrix] = {}
        self._order: List[str] = []
        for entry in d["results"]:
            lm = LabeledMatrix(entry)
            self.metrics[lm.name] = lm
            self._order.append(lm.name)
        self.n: int = d["n"]

    def __getitem__(self, name: str) -> LabeledMatrix:
        return self.metrics[name]

    def names(self) -> List[str]:
        return list(self._order)

    def __repr__(self) -> str:
        return f"EvaluationResult(metrics={self.names()}, n_extensions={self.n})"


# ---------------------------------------------------------------------------
# Linear-extension generators
# ---------------------------------------------------------------------------

class LEGenerator:
    """Exact linear-extension generator (Tree-of-Ideals enumeration)."""

    __slots__ = ("_handle", "_poset")

    def __init__(self, poset: POSet) -> None:
        self._poset = poset  # keep the source POSet alive
        self._handle = _core.build_le_generator(poset.handle)

    def get(self, from_start: bool = True, n: Optional[int] = None,
            output_every: Optional[int] = None) -> List[List[str]]:
        """Return up to ``n`` extensions (all remaining if ``n`` is None).

        Each extension is a list of element names, bottom to top.
        """
        return _core.leg_get(self._handle, bool(from_start), n, output_every)


class LEBubleyDyer:
    """Approximate linear-extension generator (Bubley–Dyer MCMC)."""

    __slots__ = ("_handle", "_poset", "seed")

    def __init__(self, poset: POSet, seed: Seed = None) -> None:
        self._poset = poset
        self._handle, self.seed = _core.build_bubley_dyer_le_generator(
            poset.handle, _seed_arg(seed))

    def get(self, from_start: bool = True, n: Optional[int] = None,
            error: Optional[float] = None,
            output_every: Optional[int] = None) -> List[List[str]]:
        """Sample linear extensions.

        Provide either ``n`` (exact number of samples) or ``error`` (target
        total-variation distance from uniform; the number of MCMC steps is
        derived). Each extension is a list of element names.
        """
        return _core.leg_bubley_dyer_get(
            self._handle, bool(from_start), n, error, output_every)


# ---------------------------------------------------------------------------
# Evaluation: MRP
# ---------------------------------------------------------------------------

def ExactMRP(poset: POSet, output_every: Optional[int] = None) -> MRPResult:
    """Exact Mutual Ranking Probability matrix (enumerates every extension)."""
    return MRPResult(_core.exact_mrp(poset.handle, output_every))


class BubleyDyerMRPGenerator:
    """Stateful MCMC estimator of the Mutual Ranking Probability matrix.

    Call :meth:`run` repeatedly to extend the same Markov chain.
    """

    __slots__ = ("_handle", "_poset", "seed")

    def __init__(self, poset: POSet, seed: Seed = None) -> None:
        self._poset = poset
        self._handle, self.seed = _core.build_bubley_dyer_mrp_generator(
            poset.handle, _seed_arg(seed))

    def run(self, n: Optional[int] = None, error: Optional[float] = None,
            output_every: Optional[int] = None) -> MRPResult:
        """Run (or continue) the chain; return the current MRP estimate."""
        return MRPResult(_core.bubley_dyer_mrp(self._handle, n, error, output_every))


# ---------------------------------------------------------------------------
# Evaluation: arbitrary metrics
# ---------------------------------------------------------------------------

def _check_metrics(functions: Sequence[str]) -> List[str]:
    names = list(functions)
    for name in names:
        if name not in VALID_METRICS:
            raise ValueError(
                f"unknown metric {name!r}; valid metrics: {', '.join(VALID_METRICS)}")
    return names


def ExactEvaluation(poset: POSet, functions: Sequence[str],
                    output_every: Optional[int] = None) -> EvaluationResult:
    """Exact evaluation of one or more metrics over every linear extension."""
    names = _check_metrics(functions)
    return EvaluationResult(_core.exact_evaluation(poset.handle, names, output_every))


class BuildBubleyDyerEvaluationGenerator:
    """Stateful MCMC evaluator of one or more metrics.

    Applies every requested metric within the same Bubley–Dyer sampling loop.
    Call :meth:`run` repeatedly to extend the same chain.
    """

    __slots__ = ("_handle", "_poset", "seed", "_functions")

    def __init__(self, poset: POSet, functions: Sequence[str], seed: Seed = None) -> None:
        self._poset = poset
        self._functions = _check_metrics(functions)
        self._handle, self.seed = _core.build_bubley_dyer_evaluation_generator(
            poset.handle, self._functions, _seed_arg(seed))

    def run(self, n: Optional[int] = None, error: Optional[float] = None,
            output_every: Optional[int] = None) -> EvaluationResult:
        """Run (or continue) the chain; return the current metric estimates."""
        return EvaluationResult(
            _core.bubley_dyer_evaluation(self._handle, n, error, output_every))


# ---------------------------------------------------------------------------
# Brüggemann–Lerche–Sørensen dominance
# ---------------------------------------------------------------------------

def BLSDominance(poset: POSet, relative: bool = False) -> Dict[str, object]:
    """Brüggemann–Lerche–Sørensen dominance matrix (absolute or relative)."""
    return _core.bls_dominance(poset.handle, bool(relative))
