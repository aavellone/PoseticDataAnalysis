"""Phase-3 full smoke test: separation, lex, fuzzy, dim-reduction, FOD."""

import poseticDataAnalysis as pda
from poseticDataAnalysis._advanced import FirstOrderDominanceAnalysis


def diamond():
    return pda.POSet(["a", "b", "c", "d"],
                     [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])


def test_separation():
    p = diamond()
    ex = pda.ExactSeparation(p, quali=["symmetric", "asymmetricLower"])
    assert set(ex.names()) == {"symmetric", "asymmetricLower"}
    gen = pda.BuildBubleyDyerSeparationGenerator(p, quali=["symmetric"], seed=3)
    res = pda.BubleyDyerSeparation(gen, n=20000)
    assert "symmetric" in res.metrics
    print("separation OK:", ex, "| MCMC n:", res.n, "seed", gen.seed)


def test_lex():
    sep = pda.LexSeparation(2, 3, types=["symmetric", "vertical"])
    assert "symmetric" in sep and "vertical" in sep and "labels" in sep
    assert len(sep["labels"]) == 9  # 3 x 3 profiles
    mrp = pda.LexMRP(2, 3)
    assert len(mrp["mrp"]) == 9
    print("lex OK: labels[:3] =", sep["labels"][:3])


def test_fuzzy():
    dom = [[1.0, 0.8, 0.3],
           [0.2, 1.0, 0.6],
           [0.7, 0.4, 1.0]]
    els = ["a", "b", "c"]
    fs = pda.FuzzySeparation(dom, els, "minimum", ["symmetric", "asymmetricLower"])
    assert "symmetric" in fs and fs["elements"] == els
    fsp = pda.FuzzySeparationProbabilistic(dom, els, ["symmetric"])
    assert "symmetric" in fsp
    fib = pda.FuzzyInBetweenness(dom, els, "minimum", ["symmetric"])
    assert len(fib["symmetric"]) == 3 and len(fib["symmetric"][0]) == 3
    print("fuzzy OK: sym[0] =", fs["symmetric"][0])


def test_dimred():
    # 4 profiles over 3 binary variables
    profile = [[0, 0, 0],
               [1, 0, 0],
               [1, 1, 0],
               [1, 1, 1]]
    weights = [1.0, 2.0, 3.0, 4.0]
    res = pda.OptimalBidimensionalEmbedding(profile, weights, "absolute")
    assert "bestLossValue" in res and "bestRepresentation" in res
    assert "variablesPriority" in res
    bpr = pda.BidimentionalPosetRepresentation(profile, weights, [1, 2, 3])
    assert "lossValue" in bpr and "representation" in bpr
    print("dimred OK: bestLoss =", res["bestLossValue"], "| bpr loss =", bpr["lossValue"])


def test_fod():
    p = pda.LinearPOSet(["1", "2", "3"])
    freq = [[10.0, 0.0],
            [5.0, 5.0],
            [0.0, 10.0]]
    res = FirstOrderDominanceAnalysis(
        p, freq, row_labels=["1", "2", "3"], col_labels=["G1", "G2"],
        metrics=["Dominance"])
    assert res["LEType"] == "TreeOfIdeals"
    dom = res["metrics"][0]
    assert dom["name"] == "Dominance"
    assert dom["cols"] == ["G1", "G2"]
    assert dom["fodClosed"] is not None
    print("FOD OK: LEType =", res["LEType"], "| le_count =", res["le_count"])
    print("       fodClosed =", dom["fodClosed"])


def test_fod_mcmc_seed():
    p = pda.POSet(["a", "b", "c", "d"], [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])
    freq = [[3.0, 1.0], [2.0, 2.0], [1.0, 3.0], [0.0, 4.0]]
    labels = ["a", "b", "c", "d"]
    r1 = FirstOrderDominanceAnalysis(p, freq, labels, ["G1", "G2"],
                                         metrics=["Dominance"], count=500, seed="42")
    r2 = FirstOrderDominanceAnalysis(p, freq, labels, ["G1", "G2"],
                                         metrics=["Dominance"], count=500, seed="42")
    assert r1["metrics"][0]["binMatrix"] == r2["metrics"][0]["binMatrix"]
    assert r1["seed"] == "42"
    print("FOD MCMC seed reproducible OK, seed =", r1["seed"])


if __name__ == "__main__":
    test_separation()
    test_lex()
    test_fuzzy()
    test_dimred()
    test_fod()
    test_fod_mcmc_seed()
    print("ALL PHASE-3 (full) TESTS PASSED")
