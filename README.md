# fmm_PB

Finite-memory fast multipole solver for the Poisson-Boltzmann equation using a Galerkin formulation. The current repository contains a CPU build and is being reorganized toward a cleaner `src/` and `include/` layout that can support later GPU work.

## Repository layout

- `src/`: solver source files
- `include/`: shared headers
- `build/`: generated object files
- `test_proteins/`: sample input cases
- `docs/`: project documentation

## Dependencies

This project should treat numerical libraries as system dependencies, not as vendored build scripts.

Required:

- C compiler (`gcc` or `clang`)
- BLAS
- LAPACK

Install notes are in [`docs/dependencies.md`](/Users/jiahuic/Garage/electrostatics/fmm_PB/docs/dependencies.md).

## Build

Configure:

```sh
cmake -S . -B build
```

Build:

```sh
cmake --build build
```

For cross-machine performance comparisons (Linux vs macOS), use a consistent
benchmark setup:

```sh
./scripts/run_apples_to_apples.sh test_proteins/1a7m
```

This script enforces:
- `Release` build
- `OpenBLAS` selection (`-DFMM_PB_BLA_VENDOR=OpenBLAS`)
- single-thread runtime (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`,
  `MKL_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`, `BLIS_NUM_THREADS=1`)

If `OpenBLAS` is unavailable, install it first (see `docs/dependencies.md`).

Clean rebuild:

```sh
rm -rf build
cmake -S . -B build
cmake --build build
```

The configure step checks BLAS and LAPACK up front and stops immediately if either is missing.

If BLAS/LAPACK live in non-default locations, pass the usual CMake search hints, for example through `CMAKE_PREFIX_PATH`.

GPU backend scaffold (CUDA, optional):

```sh
cmake -S . -B build -DFMM_PB_ENABLE_CUDA=ON
cmake --build build
```

At runtime, pass `-g=1` to request the GPU near-field backend. If the backend is unavailable or not yet fully implemented, the solver falls back to the CPU path.

## Run

```sh
./build/coulomb test_proteins/1a7m.pqr
```

## Professionalization goals

The immediate cleanup target is:

- keep generated artifacts out of git
- make external dependencies explicit and overridable
- document platform setup cleanly
- preserve a stable CPU validation path while the project structure evolves
