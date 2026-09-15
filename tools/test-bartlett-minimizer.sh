#!/usr/bin/env bash
# Run after building/installing xFitter and sourcing its dependency setup.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export LD_LIBRARY_PATH="$PWD/lib:$PWD/lib/xfitter:${LD_LIBRARY_PATH:-}"
mkdir -p temp/bartlett-minimizer-output
"${FC:-gfortran}" -cpp -ffixed-line-length-none -Iinclude -Wl,--export-dynamic \
  tools/tests/bartlett_minimizer.f -Llib -lxfitter -o temp/bartlett-minimizer-test
temp/bartlett-minimizer-test
