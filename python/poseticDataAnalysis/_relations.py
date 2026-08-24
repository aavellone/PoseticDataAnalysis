"""Binary-relation property checks and matrix-based POSet constructors.

The relation helpers mirror the R package: they operate on a raw relation given
as a set of element names plus a sequence of ``(a, b)`` pairs (meaning ``a`` is
related to ``b``), independently of any :class:`POSet`.

``POSetFromCover`` / ``POSetFromDominance`` build a :class:`POSet` from a square
0/1 matrix; they are pure-Python conveniences over :class:`POSet`, matching the
R implementation.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from . import _core
from ._poset import POSet

__all__ = [
    "IsReflexive",
    "IsSymmetric",
    "IsAntisymmetric",
    "IsTransitive",
    "IsPreorder",
    "IsPartialOrder",
    "TransitiveClosure",
    "ReflexiveClosure",
    "POSetFromCover",
    "POSetFromDominance",
]

Relation = Sequence[Tuple[str, str]]


# --- Property checks -------------------------------------------------------

def IsReflexive(set: Sequence[str], rel: Relation) -> bool:
    """True if the relation on ``set`` is reflexive (a ~ a for every a)."""
    return _core.is_reflexive(list(set), list(rel))


def IsSymmetric(rel: Relation) -> bool:
    """True if the relation is symmetric (elements inferred from the edges)."""
    return _core.is_symmetric(list(rel))


def IsAntisymmetric(rel: Relation) -> bool:
    """True if the relation is antisymmetric."""
    return _core.is_antisymmetric(list(rel))


def IsTransitive(rel: Relation) -> bool:
    """True if the relation is transitive."""
    return _core.is_transitive(list(rel))


def IsPreorder(set: Sequence[str], rel: Relation) -> bool:
    """True if the relation on ``set`` is a preorder (reflexive + transitive)."""
    return _core.is_preorder(list(set), list(rel))


def IsPartialOrder(set: Sequence[str], rel: Relation) -> bool:
    """True if the relation on ``set`` is a partial order."""
    return _core.is_partial_order(list(set), list(rel))


def TransitiveClosure(rel: Relation) -> List[Tuple[str, str]]:
    """Transitive closure of ``rel`` as a list of ``(a, b)`` pairs."""
    return _core.transitive_closure(list(rel))


def ReflexiveClosure(set: Sequence[str], rel: Relation) -> List[Tuple[str, str]]:
    """Reflexive closure of ``rel`` on ``set`` as a list of ``(a, b)`` pairs."""
    return _core.reflexive_closure(list(set), list(rel))


# --- Constructors from a square 0/1 matrix ---------------------------------

def _poset_from_square_matrix(matrix, labels, what):
    rows = [list(r) for r in matrix]
    n = len(rows)
    if any(len(r) != n for r in rows):
        raise ValueError(f"'{what}' must be a square matrix.")

    if labels is None:
        labels = [f"e{i + 1}" for i in range(n)]
    else:
        labels = list(labels)
        if len(labels) != n:
            raise ValueError(f"'labels' must have length {n}.")
        if any((lbl is None or lbl == "") for lbl in labels):
            raise ValueError("labels cannot be empty or None.")

    # Antisymmetry: no i != j with both m[i][j] and m[j][i].
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i][j] and matrix[j][i]:
                raise ValueError(
                    f"the {what} violates antisymmetry (cycle between "
                    f"{labels[i]!r} and {labels[j]!r}).")

    relations = [(labels[i], labels[j])
                 for i in range(n) for j in range(n)
                 if i != j and matrix[i][j]]
    return POSet(labels, relations)


def POSetFromCover(cover_matrix: Sequence[Sequence[int]],
                   labels: Optional[Sequence[str]] = None) -> POSet:
    """Build a POSet from a square cover (Hasse) matrix.

    ``cover_matrix[i][j]`` truthy means element ``i`` is covered by element
    ``j`` (``i < j``). ``labels`` names the elements (default ``e1, e2, ...``).
    The order relation is the transitive closure of the cover relation, computed
    by the core.
    """
    return _poset_from_square_matrix(cover_matrix, labels, "cover_matrix")


def POSetFromDominance(dominance_matrix: Sequence[Sequence[int]],
                       labels: Optional[Sequence[str]] = None) -> POSet:
    """Build a POSet from a square dominance (incidence) matrix.

    ``dominance_matrix[i][j]`` truthy means ``i <= j``. ``labels`` names the
    elements (default ``e1, e2, ...``).
    """
    return _poset_from_square_matrix(dominance_matrix, labels, "dominance_matrix")
