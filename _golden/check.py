"""Golden cross-check: compute in Python and compare against the R reference.

Reads reference/*.txt (produced by generate_reference.R) and asserts the Python
package agrees. With the same string seed the MCMC streams are bit-identical.
Run from _golden/ with the package importable (see run_compare.sh).
"""

import os
import sys

import poseticDataAnalysis as pda
from poseticDataAnalysis._advanced import FirstOrderDominanceAnalysis

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "reference")


def r_lines(name):
    with open(os.path.join(REF, name)) as f:
        return [ln.rstrip("\n") for ln in f if ln.strip() != ""]


def r_matrix(name):
    return [[float(x) for x in ln.split()] for ln in r_lines(name)]


def r_scalar(name):
    return float(r_lines(name)[0])


def maxdiff(a, b):
    return max(abs(a[i][j] - b[i][j]) for i in range(len(a)) for j in range(len(a[0])))


results = []


def check(label, ok, detail=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))


# --- Structure -------------------------------------------------------------
elems = ["a", "b", "c", "d", "e"]
dom = [("a", "c"), ("b", "c"), ("c", "e"), ("d", "e")]
p = pda.POSet(elems, dom)

check("elements", p.elements == r_lines("elements.txt"))
check("incidence matrix", [[int(v) for v in r] for r in p.incidence_matrix()]
      == [[int(v) for v in r] for r in r_matrix("incidence.txt")])
check("order relation", sorted(f"{a}<={b}" for a, b in p.order_relation()) == r_lines("order.txt"))
check("maximals", sorted(p.maximals()) == r_lines("maximal.txt"))
check("minimals", sorted(p.minimals()) == r_lines("minimal.txt"))

# --- Exact MRP -------------------------------------------------------------
emrp = pda.ExactMRP(p)
check("ExactMRP matrix", maxdiff(emrp.matrix, r_matrix("exact_mrp.txt")) < 1e-9)
check("ExactMRP n", emrp.n == int(r_lines("exact_mrp_n.txt")[0]))

# --- MCMC (bit-identical) --------------------------------------------------
py_les = pda.LEBubleyDyer(p, seed="123456789").get(from_start=True, n=25)
r_mat = r_lines("le_samples.txt")
r_les = [[r_mat[row].split()[col] for row in range(len(r_mat))]
         for col in range(len(r_mat[0].split()))]
check("MCMC linear extensions (bit-identical)", py_les == r_les,
      f"{sum(x == y for x, y in zip(py_les, r_les))}/{len(r_les)} match")

py_bmrp = pda.BubleyDyerMRPGenerator(p, seed="987654321").run(n=200000)
check("BubleyDyerMRP (bit-identical)", maxdiff(py_bmrp.matrix, r_matrix("bd_mrp.txt")) < 1e-12)

# --- Lex / Fuzzy / FOD / dim-reduction -------------------------------------
check("LexMRP(2,3)", maxdiff(pda.LexMRP(2, 3)["mrp"], r_matrix("lexmrp.txt")) < 1e-12)
check("LexSeparation(2,3) symmetric",
      maxdiff(pda.LexSeparation(2, 3, types=["symmetric"])["symmetric"], r_matrix("lexsep_sym.txt")) < 1e-12)

fdom = [[1.0, 0.2, 0.7], [0.8, 1.0, 0.4], [0.3, 0.6, 1.0]]
fs = pda.FuzzySeparationMinMax(fdom, ["a", "b", "c"], ["symmetric", "asymmetricLower"])
check("FuzzySeparation symmetric", maxdiff(fs["symmetric"], r_matrix("fuzzysep_sym.txt")) < 1e-12)
check("FuzzySeparation asymmetricLower", maxdiff(fs["asymmetricLower"], r_matrix("fuzzysep_lower.txt")) < 1e-12)

fp = pda.LinearPOSet(["1", "2", "3"])
freq = [[10.0, 0.0], [5.0, 5.0], [0.0, 10.0]]
fod = FirstOrderDominanceAnalysis(fp, freq, ["1", "2", "3"], ["G1", "G2"], metrics=["Dominance"])
dm = fod["metrics"][0]
check("FOD fodClosed", maxdiff(dm["fodClosed"], r_matrix("fod_closed.txt")) < 1e-12)
check("FOD binMatrix", maxdiff(dm["binMatrix"], r_matrix("fod_bin.txt")) < 1e-12)

profile = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1]]
dr = pda.OptimalBidimensionalEmbedding(profile, [1.0, 2.0, 3.0, 4.0], "absolute")
check("DimRed bestLossValue", abs(dr["bestLossValue"] - r_scalar("dimred_bestloss.txt")) < 1e-9)

print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
