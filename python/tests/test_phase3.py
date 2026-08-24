"""Phase-3a smoke test: relation properties, closures, matrix constructors."""

import poseticDataAnalysis as pda
from poseticDataAnalysis._relations import POSetFromCover, POSetFromDominance


def test_relation_properties():
    assert pda.IsSymmetric([("a", "b"), ("b", "a")]) is True
    assert pda.IsSymmetric([("a", "b")]) is False

    assert pda.IsAntisymmetric([("a", "b")]) is True
    assert pda.IsAntisymmetric([("a", "b"), ("b", "a")]) is False

    assert pda.IsTransitive([("a", "b"), ("b", "c"), ("a", "c")]) is True
    assert pda.IsTransitive([("a", "b"), ("b", "c")]) is False

    assert pda.IsReflexive(["a", "b"], [("a", "a"), ("b", "b")]) is True
    assert pda.IsReflexive(["a", "b"], [("a", "a")]) is False

    # Full partial order on {a,b}: a<=a, b<=b, a<=b
    assert pda.IsPartialOrder(["a", "b"], [("a", "a"), ("b", "b"), ("a", "b")]) is True
    # Not a partial order: missing reflexive pairs
    assert pda.IsPartialOrder(["a", "b"], [("a", "b")]) is False
    print("relation properties OK")


def test_closures():
    tc = set(pda.TransitiveClosure([("a", "b"), ("b", "c")]))
    assert ("a", "c") in tc and ("a", "b") in tc and ("b", "c") in tc
    assert pda.IsTransitive(list(tc)) is True

    rc = set(pda.ReflexiveClosure(["a", "b"], [("a", "b")]))
    assert ("a", "a") in rc and ("b", "b") in rc and ("a", "b") in rc
    print("closures OK: TC=", sorted(tc))


def test_poset_from_cover():
    # chain e1 < e2 < e3 as a cover matrix
    cover = [[0, 1, 0],
             [0, 0, 1],
             [0, 0, 0]]
    p = POSetFromCover(cover)
    assert p.elements == ["e1", "e2", "e3"]
    assert p.dominates("e3", "e1") == [True]      # via transitive closure
    assert p.is_comparable_with("e1", "e3") == [True]
    print("POSetFromCover OK:", p.order_relation())


def test_poset_from_dominance():
    # a <= b, a <= c ; b, c incomparable
    dom = [[1, 1, 1],
           [0, 1, 0],
           [0, 0, 1]]
    p = POSetFromDominance(dom, labels=["a", "b", "c"])
    assert p.minimals() == ["a"]
    assert p.is_incomparable_with("b", "c") == [True]
    print("POSetFromDominance OK:", p.order_relation())


def test_antisymmetry_guard():
    try:
        POSetFromDominance([[1, 1], [1, 1]])  # cycle a<=b and b<=a
    except ValueError as e:
        print("antisymmetry guard OK:", e)
    else:
        raise AssertionError("expected antisymmetry ValueError")


if __name__ == "__main__":
    test_relation_properties()
    test_closures()
    test_poset_from_cover()
    test_poset_from_dominance()
    test_antisymmetry_guard()
    print("ALL PHASE-3a TESTS PASSED")
