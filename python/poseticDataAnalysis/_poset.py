"""High-level, idiomatic Python API over the native ``_core`` extension.

Mirrors the R package: a :class:`POSet` object wraps an opaque native handle
(a PyCapsule owning the C++ ``POSetWrap``) exactly as the R S4 ``POSet`` object
wraps an external pointer. All heavy lifting happens in the C++20 core.
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Tuple, Union

from . import _core

__all__ = [
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
]

# A relation is a (from, to) pair meaning ``from < to``.
Relation = Tuple[str, str]
ElementArg = Union[str, Sequence[str]]


def _as_seq(x: ElementArg) -> List[str]:
    """Accept either a single element name or a sequence of names."""
    if isinstance(x, str):
        return [x]
    return list(x)


class POSet:
    """A partially ordered set backed by the native C++ core.

    Parameters
    ----------
    elements:
        Sequence of unique element names.
    relations:
        Optional sequence of ``(a, b)`` pairs, each meaning ``a < b``.
        The transitive/reflexive closure is computed by the core.

    A :class:`POSet` can also be created from the module-level factory
    functions (:func:`LinearPOSet`, :func:`ProductPOSet`, ...), which mirror
    the constructors of the R package.
    """

    __slots__ = ("_handle",)

    def __init__(
        self,
        elements: Optional[Sequence[str]] = None,
        relations: Optional[Sequence[Relation]] = None,
        *,
        _handle: object = None,
    ) -> None:
        if _handle is not None:
            # Internal path: wrap an existing native handle.
            self._handle = _handle
            return
        if elements is None:
            raise TypeError("POSet() requires 'elements' (a sequence of names).")
        self._handle = _core.build_poset(list(elements), list(relations or []))

    # -- internal --------------------------------------------------------
    @classmethod
    def _wrap(cls, handle: object) -> "POSet":
        return cls(_handle=handle)

    @property
    def handle(self) -> object:
        """The opaque native handle (PyCapsule). Rarely needed directly."""
        return self._handle

    # -- structure -------------------------------------------------------
    def __len__(self) -> int:
        return _core.size(self._handle)

    @property
    def elements(self) -> List[str]:
        """The element names, in internal id order."""
        return _core.elements(self._handle)

    def incidence_matrix(self) -> List[List[int]]:
        """Adjacency matrix of the order relation (rows/cols in element order)."""
        return _core.incidence_matrix(self._handle)

    def cover_matrix(self) -> List[List[int]]:
        """Cover (Hasse) matrix."""
        return _core.cover_matrix(self._handle)

    def order_relation(self) -> List[Tuple[str, str]]:
        """All strict order pairs ``(u, v)`` with ``u < v``."""
        return _core.order_relation(self._handle)

    def cover_relation(self) -> List[Tuple[str, str]]:
        """Hasse-diagram edges (covering pairs)."""
        return _core.cover_relation(self._handle)

    def incomparabilities(self) -> List[Tuple[str, str]]:
        """All incomparable element pairs."""
        return _core.incomparabilities(self._handle)

    # -- pairwise comparisons -------------------------------------------
    def dominates(self, a: ElementArg, b: ElementArg) -> List[bool]:
        """Element-wise ``a_k >= b_k``."""
        return _core.dominates(self._handle, _as_seq(a), _as_seq(b))

    def is_dominated_by(self, a: ElementArg, b: ElementArg) -> List[bool]:
        """Element-wise ``a_k <= b_k``."""
        return _core.is_dominated_by(self._handle, _as_seq(a), _as_seq(b))

    def is_comparable_with(self, a: ElementArg, b: ElementArg) -> List[bool]:
        """Element-wise comparability of ``a_k`` and ``b_k``."""
        return _core.is_comparable_with(self._handle, _as_seq(a), _as_seq(b))

    def is_incomparable_with(self, a: ElementArg, b: ElementArg) -> List[bool]:
        """Element-wise incomparability of ``a_k`` and ``b_k``."""
        return _core.is_incomparable_with(self._handle, _as_seq(a), _as_seq(b))

    # -- upset / downset -------------------------------------------------
    def upset_of(self, elements: ElementArg) -> List[str]:
        return _core.upset_of(self._handle, _as_seq(elements))

    def downset_of(self, elements: ElementArg) -> List[str]:
        return _core.downset_of(self._handle, _as_seq(elements))

    def is_upset(self, elements: ElementArg) -> bool:
        return _core.is_upset(self._handle, _as_seq(elements))

    def is_downset(self, elements: ElementArg) -> bool:
        return _core.is_downset(self._handle, _as_seq(elements))

    # -- per-element (in)comparability ----------------------------------
    def comparability_set_of(self, element: str) -> List[str]:
        return _core.comparability_set_of(self._handle, element)

    def incomparability_set_of(self, element: str) -> List[str]:
        return _core.incomparability_set_of(self._handle, element)

    # -- extremal --------------------------------------------------------
    def maximals(self) -> List[str]:
        return _core.maximals(self._handle)

    def minimals(self) -> List[str]:
        return _core.minimals(self._handle)

    def is_maximal(self, element: str) -> bool:
        return _core.is_maximal(self._handle, element)

    def is_minimal(self, element: str) -> bool:
        return _core.is_minimal(self._handle, element)

    def meet(self, elements: ElementArg) -> Optional[str]:
        """Greatest lower bound, or ``None`` if it does not exist."""
        return _core.meet(self._handle, _as_seq(elements))

    def join(self, elements: ElementArg) -> Optional[str]:
        """Least upper bound, or ``None`` if it does not exist."""
        return _core.join(self._handle, _as_seq(elements))

    # -- relations between posets ---------------------------------------
    def is_extension_of(self, other: "POSet") -> bool:
        return _core.is_extension_of(self._handle, other._handle)

    def __repr__(self) -> str:
        els = self.elements
        preview = ", ".join(els[:6]) + (", ..." if len(els) > 6 else "")
        return f"POSet(n={len(els)}, elements=[{preview}])"


# ---------------------------------------------------------------------------
# Factory functions mirroring the R constructors.
# ---------------------------------------------------------------------------

def LinearPOSet(elements: Sequence[str]) -> POSet:
    """Total order (chain) on ``elements`` in the given order."""
    return POSet._wrap(_core.build_linear_poset(list(elements)))


def _handles(posets: Iterable[POSet]) -> List[object]:
    out = []
    for p in posets:
        if not isinstance(p, POSet):
            raise TypeError("expected a sequence of POSet objects.")
        out.append(p._handle)
    return out


def ProductPOSet(posets: Iterable[POSet]) -> POSet:
    """Direct (componentwise) product order of several posets."""
    return POSet._wrap(_core.build_product_poset(_handles(posets)))


def LexicographicProductPOSet(posets: Iterable[POSet]) -> POSet:
    """Lexicographic product order (the input order is significant)."""
    return POSet._wrap(_core.build_lexicographic_product_poset(_handles(posets)))


def IntersectionPOSet(posets: Iterable[POSet]) -> POSet:
    """Intersection of the order relations (must share the same elements)."""
    return POSet._wrap(_core.build_intersection_poset(_handles(posets)))


def LinearSumPOSet(posets: Iterable[POSet]) -> POSet:
    """Linear (ordinal) sum: each poset stacked above the previous one."""
    return POSet._wrap(_core.build_linear_sum_poset(_handles(posets)))


def DisjointSumPOSet(posets: Iterable[POSet]) -> POSet:
    """Disjoint sum: the posets placed side by side, mutually incomparable."""
    return POSet._wrap(_core.build_disjoint_sum_poset(_handles(posets)))


def LiftingPOSet(poset: POSet, new_element: str) -> POSet:
    """Add ``new_element`` as a new global minimum below ``poset``."""
    return POSet._wrap(_core.build_lifting_poset(poset._handle, new_element))


def BinaryVariablePOSet(variables: Sequence[str]) -> POSet:
    """POSet of binary profiles over the given variables."""
    return POSet._wrap(_core.build_binary_variable_poset(list(variables)))


def FencePOSet(elements: Sequence[str], orientation: bool = True) -> POSet:
    """Fence (zig-zag) poset on ``elements``."""
    return POSet._wrap(_core.build_fence_poset(list(elements), bool(orientation)))


def CrownPOSet(elements_1: Sequence[str], elements_2: Sequence[str]) -> POSet:
    """Crown poset S_n^0 on two equal-length rows of elements."""
    return POSet._wrap(_core.build_crown_poset(list(elements_1), list(elements_2)))


def DualPOSet(poset: POSet) -> POSet:
    """The order-reversed (dual) poset."""
    return POSet._wrap(_core.build_dual_poset(poset._handle))
