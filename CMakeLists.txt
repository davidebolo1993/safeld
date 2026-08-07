cmake_minimum_required(VERSION 3.15)
project(safeld CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Optimising for the building CPU is off by default. On a cluster the login node
# that compiles is often not the node that runs, and a binary built with
# -march=native dies with SIGILL on any host lacking those instructions. The same
# applies to a container image built on one machine and run on another. Turn it
# on with -DSAFELD_NATIVE=ON when the build and run hosts are the same.
option(SAFELD_NATIVE "Optimise for the building CPU (faster, not portable)" OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    if(SAFELD_NATIVE)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")
        message(STATUS "SAFELD_NATIVE=ON: binary will only run on this CPU family")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Use conda environment paths
if(DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Using conda environment: $ENV{CONDA_PREFIX}")
    set(CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX};${CMAKE_PREFIX_PATH}")
    include_directories("$ENV{CONDA_PREFIX}/include")
    link_directories("$ENV{CONDA_PREFIX}/lib")
endif()

# Find packages
find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(ZLIB REQUIRED)
find_package(PkgConfig REQUIRED)

# Find htslib
pkg_check_modules(HTSLIB REQUIRED htslib)

# pkg-config reports the library directory separately from the library name.
# Without this the link works only when htslib happens to sit in a default
# search path, and fails on Homebrew/Cellar or any non-standard prefix.
if(HTSLIB_LIBRARY_DIRS)
    link_directories(${HTSLIB_LIBRARY_DIRS})
endif()

# ---------------------------------------------------------------------------
# Optional native .pgen support via plink-ng's pgenlib.
#
# pgenlib is LGPL-3.0 while safeld is MIT. Building it as a SHARED library keeps
# the LGPL relinking requirement satisfied without affecting safeld's own
# licence: a user can substitute their own pgenlib build. Do not switch this to
# a static library without reading the LGPL terms first.
#
#   git clone --depth 1 https://github.com/chrchang/plink-ng
#   cmake -DSAFELD_PGEN=ON -DPLINK_NG_DIR=/path/to/plink-ng ..
# ---------------------------------------------------------------------------
option(SAFELD_PGEN "Build native .pgen support (requires plink-ng, LGPL-3.0)" OFF)
set(PLINK_NG_DIR "" CACHE PATH "Path to a plink-ng checkout")

if(SAFELD_PGEN)
    if(NOT PLINK_NG_DIR)
        message(FATAL_ERROR "SAFELD_PGEN=ON requires -DPLINK_NG_DIR=/path/to/plink-ng")
    endif()
    set(PGENLIB_INCLUDE "${PLINK_NG_DIR}/2.0/include")
    if(NOT EXISTS "${PGENLIB_INCLUDE}/pgenlib_read.h")
        message(FATAL_ERROR "pgenlib_read.h not found under ${PGENLIB_INCLUDE}")
    endif()
    if(NOT EXISTS "${PLINK_NG_DIR}/2.0/simde")
        message(FATAL_ERROR "plink-ng's vendored simde/ is missing; clone without --sparse")
    endif()

    add_library(pgenlib SHARED
        ${PGENLIB_INCLUDE}/plink2_base.cc
        ${PGENLIB_INCLUDE}/plink2_bits.cc
        ${PGENLIB_INCLUDE}/plink2_simd.cc
        ${PGENLIB_INCLUDE}/pgenlib_misc.cc
        ${PGENLIB_INCLUDE}/pgenlib_read.cc
    )
    target_include_directories(pgenlib PUBLIC ${PGENLIB_INCLUDE})
    set_target_properties(pgenlib PROPERTIES CXX_STANDARD 11 POSITION_INDEPENDENT_CODE ON)
    message(STATUS "Native .pgen support: ON (pgenlib from ${PLINK_NG_DIR}, LGPL-3.0)")
endif()

# Shared/common source files
set(SHARED_SOURCES
    src/vcf_processor.cpp
    src/utils.cpp
)

# Main safeld executable (chunked workflow)
set(SAFELD_SOURCES
    src/main.cpp
    src/preprocessor.cpp
    src/chunked_simulator.cpp
    src/vcf_merger.cpp
    ${SHARED_SOURCES}
)

if(SAFELD_PGEN)
    list(APPEND SAFELD_SOURCES src/pgen_reader.cpp)
endif()

add_executable(safeld ${SAFELD_SOURCES})

if(SAFELD_PGEN)
    target_link_libraries(safeld pgenlib)
    target_compile_definitions(safeld PRIVATE SAFELD_HAVE_PGEN=1)
endif()

target_link_libraries(safeld 
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ZLIB::ZLIB
    OpenMP::OpenMP_CXX
    ${HTSLIB_LIBRARIES}
    pthread
    m
    dl
)

target_include_directories(safeld PRIVATE
    src/
    ${HTSLIB_INCLUDE_DIRS}
)

target_compile_options(safeld PRIVATE ${HTSLIB_CFLAGS_OTHER})

# Summary
message(STATUS "=== Build Configuration ===")
message(STATUS "HTSlib found: ${HTSLIB_FOUND}")
message(STATUS "HTSlib version: ${HTSLIB_VERSION}")
message(STATUS "Building executable:")
message(STATUS "  - safeld (three-stage chunked workflow)")
message(STATUS "Native .pgen support: ${SAFELD_PGEN}")
message(STATUS "Tune for building CPU (SAFELD_NATIVE): ${SAFELD_NATIVE}")
message(STATUS "============================")
