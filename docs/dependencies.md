# Dependencies

This project depends on external numerical libraries and should be built against installed system packages rather than ad hoc local install scripts.

## Required libraries

- BLAS
- LAPACK

The GMRES implementation is now native C and still uses BLAS/LAPACK for the underlying linear algebra kernels.

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

## Debian/Ubuntu

```sh
sudo apt-get update
sudo apt-get install build-essential libblas-dev liblapack-dev
```

Configure and build with:

```sh
cmake -S . -B build
cmake --build build
```

## Fedora/RHEL

```sh
sudo dnf install gcc blas-devel lapack-devel
```

Configure and build with:

```sh
cmake -S . -B build
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
