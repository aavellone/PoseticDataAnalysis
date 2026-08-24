# poseticDataAnalysis — unified workspace

Build and analyse **partially ordered sets (posets)** over one shared C++20
core, exposed through three self-contained front-ends that each vendor the core.

```
cpp/       Native C++20 project (Xcode + makefile). Holds the single source of
           truth for the core in cpp/PoseticDataAnalysis/POSets/ (poset engine,
           linear-extension generators, MRP / separation / dominance, fuzzy
           measures, dimensionality reduction) plus the R glue and a test main.
R/         R package (bare R C-API / .Call), vendoring the core into R/src/.
python/    Python package (bare CPython C-API — no pybind11/Cython/numpy),
           vendoring the core into python/src/core/.
_golden/   Cross-language invariance gates: the R↔Python comparison harness
           (same inputs, same 64-bit seed -> bit-identical MCMC).
sync_core.sh   Vendors the shared core into both packages (idempotent).
```

## The shared core

The R-agnostic C++20 core lives in `cpp/PoseticDataAnalysis/POSets/` (plus the
R-free handle `cpp/PoseticDataAnalysis/R/poset_wrapper.{h,cpp}`). It is **vendored**
into `R/src/` and `python/src/core/` so each package builds standalone. After
editing any core file, re-vendor both packages:

```bash
sh sync_core.sh
```

The RNG is portable **SplitMix64** (explicit bit arithmetic), so MCMC results are
bit-identical across standard libraries and platforms (Linux / Windows / macOS).

> ⚠️ Neither setuptools nor `R CMD INSTALL` rebuild a `.cpp` when only a **header**
> changes. After `sync_core.sh` (or any core-header edit) force a clean rebuild:
> `rm -rf python/build python/poseticDataAnalysis/_core*.so` and `rm -f R/src/*.o`.

## Build each front-end

- **cpp/** — `make` (see `cpp/makefile`), or open the Xcode project
  `cpp/PoseticDataAnalysis.xcodeproj`.
- **python/** — `pip install ./python` (build dep: setuptools only).
- **R/** — `R CMD INSTALL R`.

## Cross-language validation

`_golden/` runs the same computations in R and Python and checks they agree.
With the same string seed the two bindings produce **bit-identical** MCMC
streams (they share `random.h`). See `_golden/run_compare.sh`.

## Continuous integration

- `python/.github/workflows/build.yml` — builds wheels for Linux (x86_64 +
  aarch64), macOS (x86_64 + arm64) and Windows (AMD64), Python 3.9–3.13, via
  `cibuildwheel`; builds the sdist; publishes to (Test)PyPI on `v*` tags through
  Trusted Publishing (OIDC).
- `R/.github/workflows/R-CMD-check.yaml` — `R CMD check` on Linux/macOS/Windows.

## License

GPL-2.0-or-later. See `python/LICENSE`.
