"""Native build configuration for the poseticDataAnalysis Python package.

The binding is written directly against the CPython C-API (no pybind11,
no Cython, no numpy) and links the same R-agnostic C++20 core used by the
R package. Only setuptools — the standard Python build tool — is required.
"""

import glob
import os
import platform

# The C++20 core uses std::format on floating-point values, which libc++ only
# provides with a macOS 13.3+ deployment target. Set it before setuptools reads
# the environment so it is not overridden by Python's own (older) build target.
if platform.system() == "Darwin" and "MACOSX_DEPLOYMENT_TARGET" not in os.environ:
    os.environ["MACOSX_DEPLOYMENT_TARGET"] = "14.0"

from setuptools import Extension, setup

HERE = os.path.dirname(os.path.abspath(__file__))
CORE_DIR = os.path.join(HERE, "src", "core")
WRAP_DIR = os.path.join(HERE, "src", "pywrapper")

# setuptools requires '/'-separated paths RELATIVE to setup.py, never absolute.
_abs_sources = sorted(glob.glob(os.path.join(CORE_DIR, "*.cpp")))
_abs_sources += sorted(glob.glob(os.path.join(WRAP_DIR, "*.cpp")))
sources = [os.path.relpath(p, HERE).replace(os.sep, "/") for p in _abs_sources]
INCLUDE_DIRS = [
    os.path.relpath(CORE_DIR, HERE).replace(os.sep, "/"),
    os.path.relpath(WRAP_DIR, HERE).replace(os.sep, "/"),
]

# --- Platform-specific C++20 flags -----------------------------------------
if platform.system() == "Windows":
    extra_compile_args = ["/std:c++20", "/O2", "/EHsc"]
    extra_link_args = []
else:
    extra_compile_args = ["-std=c++20", "-O2"]
    extra_link_args = []
    # On macOS the deployment target is carried via MACOSX_DEPLOYMENT_TARGET
    # (set above) so the std::format float support in libc++ is available.

core_ext = Extension(
    name="poseticDataAnalysis._core",
    sources=sources,
    include_dirs=INCLUDE_DIRS,
    language="c++",
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
)

setup(ext_modules=[core_ext])
