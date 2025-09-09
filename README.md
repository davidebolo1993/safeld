# SAFELD - SAFE LD Simulator

A high-performance C++ implementation of the SAFE LD (Linkage Disequilibrium) simulator for generating synthetic genomic data from VCF files.

## Overview

SAFELD processes VCF files to generate synthetic traits while preserving the linkage disequilibrium structure of the original data. This tool is useful for:

- Generating synthetic genomic datasets for testing and validation
- Privacy-preserving genomic data sharing
- Method development and benchmarking in genomics research

## Features

- **High Performance**: Optimized C++ implementation with OpenBLAS integration
- **Memory Efficient**: Custom memory pools and optimized data structures
- **Parallel Processing**: Multi-threaded variant processing
- **Flexible Input**: Supports compressed and uncompressed VCF files
- **HTSlib Integration**: Robust VCF parsing using industry-standard library
- **Docker Support**: Containerized deployment for reproducibility

## Installation

### Prerequisites

- C++20 compatible compiler (GCC 11+ recommended)
- CMake 3.15+
- HTSlib
- OpenBLAS/LAPACK
- OpenMP

### Using Conda (Recommended)

```bash
# Create conda environment
conda env create -f safeld.yaml
conda activate safeld_conda_environment
# Build the project
mkdir build && cd build
cmake ..
make -j $(nproc)
```

### Manual Installation

```bash
sudo apt-get update
sudo apt-get install build-essential cmake libhts-dev libopenblas-dev libomp-dev
mkdir build && cd build
cmake ..
make -j $(nproc)
```

### Using Docker

```bash
# Build Docker image
docker build -t safeld .
# Run with Docker
docker run --rm -v $(pwd):/data safeld -vcf /data/input.vcf.gz -out /data/output.vcf.gz -compress
```

## Usage

### Basic Usage

```bash
# Vcf
./safeld -vcf input.vcf -out output.vcf
# Vcf.gz
./safled -vcf input.vcf.gz -out output.vcf.gz -compress
```

### Advanced Options

```txt
/.safeld [OPTIONS]

Options:
-vcf FILE Input VCF file (required, supports .vcf, .vcf.gz)
-out FILE Output VCF file (default: SAFE_LD.vcf)
-samples LIST Comma-separated sample IDs to include
-maf FLOAT Minimum allele frequency filter (default: 0.01)
-ntraits INT Number of synthetic traits (default: 10)
-workers INT Number of worker threads (default: auto-detect)
-compress Compress output (bgzip)
-h, --help Show help message
```

## Algorithm

SAFELD implements a sophisticated simulation algorithm:

1. **VCF Processing**: Parses input VCF and applies MAF filtering
2. **Duplicate Removal**: Removes duplicate variants from filtered set
3. **Traits Generation**: Creates synthetic traits matrix using random normal distribution
4. **Variant Simulation**: 
   - Standardizes original dosages per variant
   - Computes synthetic dosages via matrix multiplication
   - Scales results to valid dosage range [0, 2]
5. **Output Generation**: Writes synthetic VCF with preserved structure

## Performance

Benchmarks on chromosome 22 data (subset of UKBiobank genotype data):
- **Input**: 134,987 variants, 4,783 samples
- **After filtering**: 134,987 variants, 4,783 samples
- **Processing time**: 17m 11s   
- **Memory usage**: 16,05 Gb
- **Throughput**: XXX

## Technical Details

### Architecture

- **VCFProcessor**: HTSlib-based VCF parsing and filtering
- **SimulationEngine**: Matrix operations and synthetic data generation
- **MemoryPool**: Optimized memory management for high-throughput processing
- **Utils**: Logging, timing, and mathematical utilities

### Dependencies

- **HTSlib**: VCF/BCF file format handling
- **OpenBLAS**: Optimized linear algebra operations
- **OpenMP**: Parallel processing support
- **BGZF**: Block gzip compression for output

## File Formats

### Input VCF Requirements

- Must contain `DS` (dosage) format field
- Should include `AF` (allele frequency) in INFO field (calculated if missing)
- Supports both compressed (.vcf.gz) and uncompressed (.vcf) files

### Output VCF Structure

```txt
##fileformat=VCFv4.1
##source=safeld-cpp
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT T1 T2 ...
chr1 1000 rs123 A G . PASS . DS 1.23 0.45 ..
```
## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use SAFELD in your research, please cite:
XXX
