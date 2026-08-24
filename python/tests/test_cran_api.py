"""Validate the CRAN-mirroring public function API."""

import poseticDataAnalysis as pda


def diamond():
    return pda.POSet(["a", "b", "c", "d"],
                     [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])


def test_free_query_functions():
    p = diamond()
    assert pda.POSetElements(p) == ["a", "b", "c", "d"]
    assert pda.POSetMinimals(p) == ["a"]
    assert pda.POSetMaximals(p) == ["d"]
    assert pda.Dominates(p, "d", "a") == [True]
    assert pda.IsDominatedBy(p, "a", "d") == [True]
    assert pda.IsComparableWith(p, "b", "c") == [False]
    assert pda.IsIncomparableWith(p, "b", "c") == [True]
    assert set(pda.UpsetOf(p, "a")) == {"a", "b", "c", "d"}
    assert set(pda.DownsetOf(p, "d")) == {"a", "b", "c", "d"}
    assert pda.IsUpset(p, ["d"]) is True
    assert pda.POSetMeet(p, ["b", "c"]) == "a"
    assert pda.POSetJoin(p, ["b", "c"]) == "d"
    assert pda.IsMaximal(p, "d") and pda.IsMinimal(p, "a")
    assert len(pda.DominanceMatrix(p)) == 4
    assert pda.OrderRelation(p)
    assert pda.CoverRelation(p)
    assert pda.IncomparabilityRelation(p) == [("b", "c")] or ("b", "c") in pda.IncomparabilityRelation(p)
    assert pda.ComparabilitySetOf(p, "a")
    assert "d" in pda.IncomparabilitySetOf(p, "b") or pda.IncomparabilitySetOf(p, "b") == ["c"]
    assert pda.IsExtensionOf(pda.LinearPOSet(["a", "b", "c", "d"]), p) is True
    print("free query functions OK")


def test_le_and_evaluation_drivers():
    p = diamond()
    # exact
    les = pda.LEGet(pda.LEGenerator(p))
    assert len(les) == 2
    # MCMC
    s = pda.LEGet(pda.LEBubleyDyer(p, seed="7"), n=50)
    assert len(s) == 50
    r = pda.BubleyDyerMRP(pda.BubleyDyerMRPGenerator(p, seed="7"), n=10000)
    assert r.n >= 10000
    e = pda.BubleyDyerEvaluation(
        pda.BuildBubleyDyerEvaluationGenerator(p, ["MutualRankingProbability"], seed="7"),
        n=10000)
    assert "MutualRankingProbability" in e.metrics
    print("LE / evaluation drivers OK")


def test_public_surface_matches_cran():
    # The non-CRAN functions must NOT be part of the public API.
    for name in ("FirstOrderDominanceAnalysis", "POSetFromCover", "POSetFromDominance"):
        assert not hasattr(pda, name), f"{name} should be hidden from the public API"
    # A representative sample of the CRAN functions must be present.
    for name in ("POSet", "OrderRelation", "Dominates", "ExactMRP", "BubleyDyerMRP",
                 "LEGet", "FuzzySeparationMinMax", "OptimalBidimensionalEmbedding",
                 "TransitiveClosure", "POSetMeet"):
        assert hasattr(pda, name), f"{name} must be public"
    print(f"public API surface OK ({len(pda.__all__)} names in __all__)")


if __name__ == "__main__":
    test_free_query_functions()
    test_le_and_evaluation_drivers()
    test_public_surface_matches_cran()
    print("ALL CRAN-API TESTS PASSED")
