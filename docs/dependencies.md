# Dependencies

This project depends on external numerical libraries and should be built against installed system packages rather than ad hoc local install scripts.

## Required libraries

- BLAS
- LAPACK

## Optional library

- `libf2c`

`libf2c` is only needed on systems where the translated GMRES routine requires it at link time. The repository includes a compatibility header at [`include/f2c.h`](/Users/jiahuic/Garage/electrostatics/fmm_PB/include/f2c.h), but the runtime/library dependency should come from your package manager if needed.

## macOS

If the default Accelerate-backed BLAS/LAPACK setup works, a plain build is usually enough:

```sh
make
```

If you prefer Homebrew OpenBLAS/LAPACK, install packages first and then override build variables as needed:

```sh
brew install openblas lapack
make BLAS_LIBS=-lopenblas LAPACK_LIBS=-llapack
```

If `libf2c` is required on your machine, install it separately and pass:

```sh
make F2C_LIBS=-lf2c
```

## Debian/Ubuntu

```sh
sudo apt-get update
sudo apt-get install build-essential libblas-dev liblapack-dev libf2c2-dev
```

Configure and build with:

```sh
cmake -S . -B build -DFMM_PB_USE_F2C=ON
cmake --build build
```

## Fedora/RHEL

```sh
sudo dnf install gcc blas-devel lapack-devel libf2c
```

Configure and build with:

```sh
cmake -S . -B build -DFMM_PB_USE_F2C=ON
cmake --build build
```

## Verifying linkage

Inspect detected libraries during CMake configure:

```sh
cmake -S . -B build
```

Then build:

```sh
cmake --build build
```
