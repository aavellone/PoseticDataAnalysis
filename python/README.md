# poseticDataAnalysis (Python)

Build and analyse **partially ordered sets (posets)** in Python, backed by the
same R-agnostic **C++20 core** used by the R package of the same name.

The binding is written **directly against the CPython C-API** — no pybind11,
no Cython, no numpy, no third-party runtime dependencies. It mirrors the R
package architecture: a thin native layer over a shared HPC core.

## Architecture

```
src/core/       R-agnostic C++20 core (shared with the R package)
src/pywrapper/  CPython C-API binding (PyObject* <-> C++), the analog of the R rwrapper_*
poseticDataAnalysis/  idiomatic Python layer (the POSet class + constructors)
```

A `POSet` object holds an opaque native handle (a `PyCapsule` owning the C++
`POSetWrap`), exactly as the R S4 `POSet` holds an external pointer. The capsule
frees the C++ object when Python's garbage collector reclaims it.

## Build & install

Requires a C++20 compiler and CPython ≥ 3.9 (macOS deployment target 14.0).

```bash
# in-place build for development
python3 setup.py build_ext --inplace

# or a normal install
python3 -m pip install .
```

## Quick start

```python
import poseticDataAnalysis as pda

# Diamond poset: a < b, a < c, b < d, c < d
p = pda.POSet(["a", "b", "c", "d"],
              [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])

len(p)                       # 4
p.elements                   # ['a', 'b', 'c', 'd']
p.minimals(), p.maximals()   # (['a'], ['d'])
p.is_comparable_with("b", "c")   # [False]
p.upset_of("a")                  # ['a', 'b', 'c', 'd']
p.meet(["b", "c"]), p.join(["b", "c"])   # ('a', 'd')

# Constructors mirror the R package
chain = pda.LinearPOSet(["x", "y", "z"])
prod  = pda.ProductPOSet([pda.LinearPOSet(["a", "b"]),
                          pda.LinearPOSet(["1", "2"])])
dual  = pda.DualPOSet(chain)
```

## Public API

The published package exposes **exactly the functions of the CRAN R package**
`poseticDataAnalysis` (same names, functional form: `Dominates(poset, a, b)`,
`OrderRelation(poset)`, `POSetElements(poset)`, …). `poseticDataAnalysis.__all__`
matches the CRAN NAMESPACE one-to-one (66 functions).

Implemented but intentionally **not** part of the public API for now (reachable
only via submodules, e.g. `poseticDataAnalysis._advanced`):
`FirstOrderDominanceAnalysis`, `POSetFromCover`, `POSetFromDominance`. Method-style
access on the `POSet` object (`poset.dominates(...)`) also exists as a convenience.

## Status

Implemented and tested:

- POSet constructors: `POSet`, `LinearPOSet`, `ProductPOSet`,
  `LexicographicProductPOSet`, `IntersectionPOSet`, `LinearSumPOSet`,
  `DisjointSumPOSet`, `LiftingPOSet`, `BinaryVariablePOSet`, `FencePOSet`,
  `CrownPOSet`, `DualPOSet`.
- Structure & queries: elements, incidence/cover matrices, order/cover
  relations, incomparabilities, dominance/comparability (vectorized), up/down
  sets, comparability sets, extremal elements, meet/join, extension check.
- Linear-extension generators: `LEGenerator` (exact, Tree-of-Ideals) and
  `LEBubleyDyer` (MCMC), 64-bit seeds carried as decimal strings.
- Evaluation: `ExactMRP` / `BubleyDyerMRPGenerator` (Mutual Ranking
  Probability, exact and MCMC), `ExactEvaluation` /
  `BuildBubleyDyerEvaluationGenerator` (multiple metrics:
  MutualRankingProbability, AverageHeight, symmetric, asymmetricLower,
  asymmetricUpper), and `BLSDominance` (Brüggemann–Lerche–Sørensen).

```python
import poseticDataAnalysis as pda
p = pda.POSet(["a", "b", "c", "d"],
              [("a", "b"), ("a", "c"), ("b", "d"), ("c", "d")])

pda.LEGenerator(p).get()          # exact linear extensions
g = pda.LEBubleyDyer(p, seed="42")
g.get(n=1000); g.seed             # MCMC samples; seed echoed back as a string
pda.ExactMRP(p).matrix            # exact mutual ranking probability
pda.BubleyDyerMRPGenerator(p, seed=1).run(n=100000).matrix   # MCMC estimate
```

- Relation helpers: `IsReflexive`, `IsSymmetric`, `IsAntisymmetric`,
  `IsTransitive`, `IsPreorder`, `IsPartialOrder`, `TransitiveClosure`,
  `ReflexiveClosure`; matrix constructors `POSetFromCover`,
  `POSetFromDominance`.
- Separation: `ExactSeparation`, `BuildBubleyDyerSeparationGenerator` /
  `BubleyDyerSeparation`, lexicographic `LexSeparation` / `LexMRP`.
- Fuzzy measures on a dominance matrix: `FuzzySeparation`,
  `FuzzyInBetweenness` (+ `MinMax` / `Probabilistic` variants).
- Dimensionality reduction: `OptimalBidimensionalEmbedding`,
  `BidimentionalPosetRepresentation`.
- `FirstOrderDominanceAnalysis` (single/multi-poset, exact / MCMC / multi-chain).

Not exposed in the Python binding: metrics supplied as user callables (the R
binding accepts R-function metrics and R-function fuzzy norms; the Python
binding exposes the built-in metrics/norms only), and the derived fuzzy
`vertical` / `horizontal` separation types.

### Cross-validation against the R package

Every implemented function is validated against the R package on shared inputs;
with the same 64-bit string seed the two bindings produce **bit-identical** MCMC
results (they share the same C++ core and `random.h`).

### Build hygiene when the core changes

Neither setuptools nor `R CMD INSTALL` rebuild a `.cpp` when only a **header**
it includes changes. After editing a core header (e.g. `random.h`), force a full
recompile so stale objects do not linger:

```bash
# Python
rm -rf build poseticDataAnalysis/_core*.so && python3 setup.py build_ext --inplace
# R package
rm -f src/*.o src/*.so && R CMD INSTALL .
```

## Keeping the core in sync

`src/core/` is a copy of the shared C++ core (`POSets/` + `poset_wrapper.{h,cpp}`
from the Xcode project). When the core changes, re-copy those files — the same
manual-sync workflow the R package uses.
