#!/bin/sh
# Cross-language invariance gate: generate the R reference, then check that the
# Python package agrees (bit-identical MCMC with the same string seed).
#
# Requires: the R package installed (R CMD INSTALL ../R) and the Python
# extension built (cd ../python && python3 setup.py build_ext --inplace).
set -eu

cd "$(dirname "$0")"
PY_PKG="$(cd ../python && pwd)"

echo "==> Generating R reference..."
Rscript generate_reference.R

echo "==> Checking Python against reference..."
PYTHONPATH="$PY_PKG" python3 check.py
