"""End-to-end smoke test for the native poseticDataAnalysis binding."""

import poseticDataAnalysis as pda


def test_basic_poset():
    # Diamond: a < b, a < c, b < d, c < d
    p = pda.POSet(["a", "b", "c", "d"], [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])
    assert len(p) == 4
    assert set(p.elements) == {"a", "b", "c", "d"}
    assert p.minimals() == ["a"]
    assert p.maximals() == ["d"]
    assert p.dominates("d", "a") == [True]
    assert p.is_dominated_by("a", "d") == [True]
    assert p.is_comparable_with("b", "c") == [False]
    assert p.is_incomparable_with("b", "c") == [True]
    assert set(p.upset_of("a")) == {"a", "b", "c", "d"}
    assert set(p.downset_of("d")) == {"a", "b", "c", "d"}
    assert p.meet(["b", "c"]) == "a"
    assert p.join(["b", "c"]) == "d"
    assert p.is_maximal("d") and p.is_minimal("a")
    print("diamond:", repr(p))


def test_pairwise_vectorized():
    p = pda.POSet(["a", "b", "c", "d"], [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])
    res = p.is_comparable_with(["a", "b", "a"], ["d", "c", "b"])
    assert res == [True, False, True]


def test_constructors():
    chain = pda.LinearPOSet(["x", "y", "z"])
    assert chain.order_relation()  # non-empty
    assert chain.is_comparable_with("x", "z") == [True]

    prod = pda.ProductPOSet([pda.LinearPOSet(["a", "b"]), pda.LinearPOSet(["1", "2"])])
    assert len(prod) == 4

    dual = pda.DualPOSet(chain)
    assert dual.dominates("x", "z") == [True]  # order reversed


def test_error_forwarding():
    p = pda.POSet(["a", "b"], [("a", "b")])
    try:
        p.is_maximal("nope")
    except Exception as e:  # core raises -> ValueError/RuntimeError
        print("expected error:", type(e).__name__, e)
    else:
        raise AssertionError("expected an error for unknown element")


if __name__ == "__main__":
    test_basic_poset()
    test_pairwise_vectorized()
    test_constructors()
    test_error_forwarding()
    print("ALL SMOKE TESTS PASSED")
