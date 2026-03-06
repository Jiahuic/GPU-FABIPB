#!/usr/bin/env sh
set -eu

BUILD_DIR="${BUILD_DIR:-build-a2a}"

cmake -S . -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DFMM_PB_BLA_VENDOR=OpenBLAS

cmake --build "$BUILD_DIR"

# Pin common BLAS/OpenMP runtimes to one thread for reproducible comparisons.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <panel-base-or-pqr-path> [solver options...]" >&2
  echo "Example: $0 test_proteins/1a7m" >&2
  exit 2
fi

panel="$1"
shift
case "$panel" in
  *.pqr) panel="${panel%.pqr}" ;;
esac

exec "$BUILD_DIR/coulomb" "$panel" "$@"
