cmake_minimum_required(VERSION 3.16)

project(fmm_PB VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(FMM_PB_USE_F2C "Link against libf2c explicitly" OFF)

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

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

target_link_libraries(coulomb
  PRIVATE
    LAPACK::LAPACK
    BLAS::BLAS
)

if(UNIX AND NOT APPLE)
  target_link_libraries(coulomb PRIVATE m)
endif()

if(FMM_PB_USE_F2C)
  find_library(F2C_LIBRARY NAMES f2c REQUIRED)
  target_link_libraries(coulomb PRIVATE ${F2C_LIBRARY})
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(coulomb PRIVATE -g)
endif()

message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
message(STATUS "Use libf2c: ${FMM_PB_USE_F2C}")
