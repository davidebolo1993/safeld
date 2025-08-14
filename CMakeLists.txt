cmake_minimum_required(VERSION 3.15)
project(safeld CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler optimizations  
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -DNDEBUG")
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

# Core source files
set(SOURCES
    src/main.cpp
    src/vcf_processor.cpp
    src/simulation_engine.cpp
    src/memory_pool.cpp
    src/utils.cpp
)

# Create executable
add_executable(safeld ${SOURCES})

# Link libraries
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

# Include directories
target_include_directories(safeld PRIVATE 
    src/
    ${HTSLIB_INCLUDE_DIRS}
)

# Compiler flags
target_compile_options(safeld PRIVATE ${HTSLIB_CFLAGS_OTHER})

message(STATUS "=== Build Configuration ===")
message(STATUS "HTSlib found: ${HTSLIB_FOUND}")
message(STATUS "HTSlib version: ${HTSLIB_VERSION}")
message(STATUS "============================")

