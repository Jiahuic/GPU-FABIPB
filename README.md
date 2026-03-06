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

Optional, platform-dependent:

- `libf2c` for the translated GMRES routine on systems where BLAS/LAPACK packages do not already provide the required symbols

The repository currently keeps a local copy of [`include/f2c.h`](/Users/jiahuic/Garage/electrostatics/fmm_PB/include/f2c.h) as a compatibility header for the translated source in [`src/gmres.c`](/Users/jiahuic/Garage/electrostatics/fmm_PB/src/gmres.c). That is only a header, not a bundled install script or vendored binary library.

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

Clean rebuild:

```sh
rm -rf build
cmake -S . -B build
cmake --build build
```

If your system requires `libf2c`, require it during configure:

```sh
cmake -S . -B build -DFMM_PB_REQUIRE_F2C=ON
cmake --build build
```

The configure step now checks BLAS and LAPACK up front and stops immediately if either is missing. `libf2c` is detected automatically when present, and can be made mandatory with `FMM_PB_REQUIRE_F2C=ON`.

If BLAS/LAPACK live in non-default locations, pass the usual CMake search hints, for example through `CMAKE_PREFIX_PATH`.

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
