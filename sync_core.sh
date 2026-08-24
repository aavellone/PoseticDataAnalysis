#!/bin/sh
# Vendor the shared, R-agnostic C++20 core into the R and Python packages so each
# stays self-contained (see README, "The shared core ...").
#
# Single source of truth:
#   cpp/PoseticDataAnalysis/POSets/            -> the whole R-free core (*.h, *.cpp)
#   cpp/PoseticDataAnalysis/R/poset_wrapper.*  -> the R-free handle wrapper
#
# Destinations (core only; the language-specific wrappers and build artifacts
# are left untouched):
#   R/src/            (alongside rwrapper_*.cpp, r_display.*, ... which are R-only)
#   python/src/core/  (alongside src/pywrapper/ which is Python-only)
#
# Idempotent: re-running with no upstream change overwrites identical files.
# Run after editing any core header/source, then rebuild the packages
# (remember: neither setuptools nor R CMD INSTALL rebuild a .cpp on header-only
# changes, so force a clean rebuild — see README).
set -eu

ROOT="$(cd "$(dirname "$0")" && pwd)"
CORE="$ROOT/cpp/PoseticDataAnalysis/POSets"
WRAP="$ROOT/cpp/PoseticDataAnalysis/R"

for DST in "$ROOT/R/src" "$ROOT/python/src/core"; do
    mkdir -p "$DST"
    cp "$CORE"/*.h "$CORE"/*.cpp "$DST"/
    cp "$WRAP"/poset_wrapper.h "$WRAP"/poset_wrapper.cpp "$DST"/
    echo "Synced core -> ${DST#"$ROOT/"}"
done
