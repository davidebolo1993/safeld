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

add_executable(safeld ${SAFELD_SOURCES})

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
message(STATUS "============================")
