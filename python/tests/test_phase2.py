"""Phase-2 smoke test: generators + evaluation."""

import poseticDataAnalysis as pda


def diamond():
    return pda.POSet(["a", "b", "c", "d"],
                     [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])


def test_exact_le_generator():
    p = diamond()
    gen = pda.LEGenerator(p)
    les = gen.get()
    # Diamond has exactly 2 linear extensions: a,b,c,d and a,c,b,d
    assert len(les) == 2, les
    for le in les:
        assert le[0] == "a" and le[-1] == "d"
    print("exact LE:", les)


def test_exact_mrp():
    p = diamond()
    res = pda.ExactMRP(p)
    assert res.elements == ["a", "b", "c", "d"]
    assert len(res.matrix) == 4 and len(res.matrix[0]) == 4
    assert res.n == 2  # 2 linear extensions
    # a is below everything: P(a <= x) == 1 for all x
    idx = {e: i for i, e in enumerate(res.elements)}
    a = idx["a"]
    for x in range(4):
        assert abs(res.matrix[a][x] - 1.0) < 1e-9, res.matrix[a]
    print("exact MRP:", res)


def test_bubley_dyer_seed_roundtrip():
    p = diamond()
    g1 = pda.LEBubleyDyer(p, seed="424242")
    s1 = g1.get(n=6)
    assert g1.seed == "424242"
    g2 = pda.LEBubleyDyer(p, seed="424242")
    s2 = g2.get(n=6)
    assert s1 == s2, "same seed must reproduce the same samples"
    # int seed accepted too
    g3 = pda.LEBubleyDyer(p, seed=424242)
    assert g3.get(n=6) == s1
    print("seed roundtrip OK, seed =", g1.seed)


def test_bubley_dyer_mrp_converges():
    p = diamond()
    exact = pda.ExactMRP(p)
    mc = pda.BubleyDyerMRPGenerator(p, seed=1)
    approx = mc.run(n=200000)
    # Compare matrices: should be close for a 4-element poset
    max_err = max(abs(exact.matrix[i][j] - approx.matrix[i][j])
                  for i in range(4) for j in range(4))
    assert max_err < 0.02, f"MCMC MRP too far from exact: {max_err}"
    print(f"MCMC MRP max abs error vs exact: {max_err:.4f} (seed={mc.seed})")


def test_exact_evaluation():
    p = diamond()
    res = pda.ExactEvaluation(p, ["MutualRankingProbability", "AverageHeight"])
    assert set(res.names()) == {"MutualRankingProbability", "AverageHeight"}
    mrp = res["MutualRankingProbability"]
    assert len(mrp.matrix) == 4
    print("exact evaluation:", res, "| AverageHeight cols:", res["AverageHeight"].cols)


def test_bubley_dyer_evaluation():
    p = diamond()
    gen = pda.BuildBubleyDyerEvaluationGenerator(
        p, ["MutualRankingProbability"], seed=7)
    res = gen.run(n=50000)
    assert "MutualRankingProbability" in res.metrics
    assert res.n >= 50000
    print("MCMC evaluation:", res, "| seed:", gen.seed)


def test_bls_dominance():
    p = diamond()
    abs_d = pda.BLSDominance(p, relative=False)
    rel_d = pda.BLSDominance(p, relative=True)
    assert len(abs_d["matrix"]) == 4
    assert rel_d["elements"] == ["a", "b", "c", "d"]
    print("BLS absolute[0]:", abs_d["matrix"][0])


if __name__ == "__main__":
    test_exact_le_generator()
    test_exact_mrp()
    test_bubley_dyer_seed_roundtrip()
    test_bubley_dyer_mrp_converges()
    test_exact_evaluation()
    test_bubley_dyer_evaluation()
    test_bls_dominance()
    print("ALL PHASE-2 TESTS PASSED")
