cmake_minimum_required(VERSION 3.16)

project(fmm_PB VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(FMM_PB_REQUIRE_F2C "Fail configure unless libf2c is available" OFF)

set(FMM_PB_SOURCES
  src/coulomb.c
  src/gkSetup.c
  src/input.c
  src/gkGlobal.c
  src/expan.c
  src/moments.c
  src/numQuad.c
  src/fmm.c
  src/kernel.c
  src/precond_fmm.c
  src/gmres.c
  src/treecode.c
)

add_executable(coulomb ${FMM_PB_SOURCES})

target_include_directories(coulomb
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

find_package(BLAS QUIET)
if(NOT BLAS_FOUND)
  message(FATAL_ERROR
    "BLAS was not found. This project requires an installed BLAS library. "
    "libf2c is not a replacement for BLAS because src/gmres.c calls BLAS "
    "routines directly."
  )
endif()

find_package(LAPACK QUIET)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR
    "LAPACK was not found. Please install LAPACK and rerun CMake."
  )
endif()

target_link_libraries(coulomb
  PRIVATE
    LAPACK::LAPACK
    BLAS::BLAS
)

if(UNIX AND NOT APPLE)
  target_link_libraries(coulomb PRIVATE m)
endif()

find_library(F2C_LIBRARY NAMES f2c)

if(F2C_LIBRARY)
  target_link_libraries(coulomb PRIVATE ${F2C_LIBRARY})
elseif(FMM_PB_REQUIRE_F2C)
  message(FATAL_ERROR
    "FMM_PB_REQUIRE_F2C is ON, but libf2c was not found."
  )
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(coulomb PRIVATE -g)
endif()

message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
message(STATUS "libf2c library: ${F2C_LIBRARY}")
message(STATUS "Require libf2c: ${FMM_PB_REQUIRE_F2C}")
