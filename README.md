# SAFELD - SAFE LD Simulator

A high-performance C++ implementation of the SAFE LD (Linkage Disequilibrium) simulator for generating synthetic genomic data from VCF files.

## Overview

SAFELD processes VCF files to generate synthetic traits while preserving the linkage disequilibrium structure of the original data. This tool is useful for:

- Generating synthetic genomic datasets for testing and validation
- Privacy-preserving genomic data sharing
- Method development and benchmarking in genomics research

## Features

- **High Performance**: Optimized C++ implementation with OpenBLAS integration and BLAS GEMM operations
- **Memory Efficient**: Custom memory pools and optimized data structures
- **Parallel Processing**: Multi-threaded variant processing with OpenMP
- **Flexible Input**: Supports compressed and uncompressed VCF files
- **HTSlib Integration**: Robust VCF parsing using industry-standard library
- **Docker Support**: Containerized deployment for reproducibility
- **Two Operational Modes**: 
  - Standard mode for typical datasets
  - Chunked mode for very large datasets (>100K samples or high trait counts)

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

This builds two executables:
- `safeld` - Standard mode (fast, in-memory processing)
- `safeld_chunked` - Chunked mode (memory-efficient, scalable to very large datasets)

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

### Standard Mode (Recommended for most use cases)

```bash
# Basic usage with VCF
./safeld -vcf input.vcf -out output.vcf

# Compressed VCF with custom parameters
./safeld -vcf input.vcf.gz -out output.vcf.gz -compress -ntraits 5000 -maf 0.01
```

**Standard Mode Options:**

```txt
./safeld [OPTIONS]

Options:
  -vcf FILE        Input VCF file (required, supports .vcf, .vcf.gz)
  -out FILE        Output VCF file (default: SAFE_LD.vcf)
  -samples LIST    Comma-separated sample IDs to include
  -maf FLOAT       Minimum allele frequency filter (default: 0.01)
  -ntraits INT     Number of synthetic traits (default: 10)
  -workers INT     Number of worker threads (default: auto-detect)
  -compress        Compress output (bgzip)
  -h, --help       Show help message
```

### Chunked Mode (For very large datasets)

The chunked mode uses a **three-stage workflow**:

#### Stage 1: Preprocessing

Generates trait matrix and partitions variants into chunks.

```bash
./safeld_chunked preprocess \
  -vcf input.vcf.gz \
  -out preprocessed_data \
  -ntraits 10000 \
  -chunk-size 5000 \
  -maf 0.01
```

**Preprocessing Options:**

```txt
./safeld_chunked preprocess [OPTIONS]

Options:
  -vcf FILE            Input VCF file (required)
  -out DIR             Output directory for preprocessed data
  -samples LIST        Comma-separated sample IDs
  -maf FLOAT           MAF filter (default: 0.01)
  -ntraits INT         Number of traits (default: 10)
  -chunk-size INT      Variants per chunk (default: 10000)
  -traits-per-tile INT Traits per tile (default: auto, ~1GB tiles)
  -h, --help           Show this help message
```

**Output structure:**
```
preprocessed_data/
├── traits/
│   ├── W_tile_0.bin      # Trait matrix (tiled if large)
│   ├── W_tile_1.bin
│   └── metadata.txt      # Trait dimensions and tile info
└── chunks/
    ├── chunk_0.bin       # Standardized genotypes
    ├── chunk_0.meta      # Variant metadata
    ├── chunk_1.bin
    ├── chunk_1.meta
    └── ...
```

#### Stage 2: Simulation

Processes chunks to generate synthetic traits. Each chunk uses all available cores via optimized BLAS GEMM operations.

```bash
# Process all chunks sequentially
./safeld_chunked simulate \
  -prep preprocessed_data \
  -out results \
  -workers 32 \
  -compress

# Process specific chunk range (useful for cluster parallelization)
./safeld_chunked simulate \
  -prep preprocessed_data \
  -out results \
  -start-chunk 0 \
  -end-chunk 9 \
  -workers 32 \
  -compress
```

**Simulation Options:**

```txt
./safeld_chunked simulate [OPTIONS]

Options:
  -prep DIR          Preprocessed data directory (required)
  -out DIR           Output directory for results (required)
  -workers INT       Number of threads (default: auto-detect)
  -compress          Compress output VCF chunks
  -start-chunk INT   First chunk to process (default: all)
  -end-chunk INT     Last chunk to process (default: all)
  -h, --help         Show this help message
```


#### Stage 3: Merge

Combines all chunk VCF files into a single output file.

```bash
./safeld_chunked merge \
  -in results \
  -out final_output.vcf.gz
```

**Merge Options:**

```txt
./safeld_chunked merge [OPTIONS]

Options:
  -in DIR            Directory with chunk VCF files (required)
  -out FILE          Output merged VCF file (required)
  -no-compress       Don't compress output (default: compressed)
  -h, --help         Show this help message
```

### Complete Chunked Workflow Example

```bash
# 1. Preprocess (one time per dataset)
./safeld_chunked preprocess \
  -vcf large_dataset.vcf.gz \
  -out preprocessed_chr22 \
  -ntraits 10000 \
  -chunk-size 5000 \
  -maf 0.01

# 2. Simulate (can be parallelized)
./safeld_chunked simulate \
  -prep preprocessed_chr22 \
  -out results_chr22 \
  -workers 32 \
  -compress

# 3. Merge
./safeld_chunked merge \
  -in results_chr22 \
  -out chr22_final.vcf.gz
```

## Algorithm

### Standard Mode Algorithm

SAFELD implements a sophisticated simulation algorithm:

1. **VCF Processing**: Parses input VCF and applies MAF filtering
2. **Duplicate Removal**: Removes duplicate variants from filtered set
3. **Traits Generation**: Creates synthetic traits matrix using random normal distribution
4. **Variant Simulation**: 
   - Standardizes original dosages per variant
   - Computes synthetic dosages via matrix multiplication with BLAS operations
   - Scales results to valid dosage range [0, 2]
5. **Output Generation**: Writes synthetic VCF with preserved structure

### Chunked Mode Algorithm

The chunked mode separates preprocessing from simulation for scalability:

**Preprocessing Stage:**
1. Parse VCF once and apply MAF filtering
2. Generate trait matrix W (T × S) with standard normal random values
3. Standardize genotype dosages and partition into chunks (B variants each)
4. Serialize traits (tiled if large) and genotype chunks to disk

**Simulation Stage:**
1. Load trait matrix W into memory once
2. For each chunk:
   - Load standardized genotypes G (B × S)
   - Compute Y = G × W^T using optimized BLAS GEMM
   - Scale synthetic dosages to [0, 2] range
   - Write chunk VCF to disk
3. Release chunk memory before processing next

**Merge Stage:**
- Concatenate all chunk VCF files maintaining chromosome order
- Stream processing for memory efficiency


### Dependencies

- **HTSlib**: VCF/BCF file format handling
- **OpenBLAS**: Optimized linear algebra operations (cblas_ddot, cblas_dgemm)
- **OpenMP**: Parallel processing support
- **BGZF**: Block gzip compression for output

## File Formats

### Input VCF Requirements

- Must contain `DS` (dosage) format field or `GT` (genotype) field
- Should include `AF` (allele frequency) in INFO field (calculated if missing)
- Supports both compressed (.vcf.gz) and uncompressed (.vcf) files

### Output VCF Structure

```txt
##fileformat=VCFv4.1
##source=safeld-cpp
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT  T1    T2    ...
chr1    1000    rs123   A    G    .     PASS    .     DS      1.23  0.45  ...
```

### Chunked Mode Binary Formats

**Trait Matrix Tiles** (`W_tile_*.bin`):
- Row-major double-precision matrix
- Each file contains a subset of traits × all samples
- Dimensions stored in `metadata.txt`

**Genotype Chunks** (`chunk_*.bin`):
- Row-major double-precision matrix
- Standardized dosages: (dosage - mean) / std
- Dimensions: chunk_size × n_samples

**Chunk Metadata** (`chunk_*.meta`):
- Text format with variant annotations
- Format: CHROM POS ID REF ALT (one variant per line)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use SAFELD in your research, please cite:
[Citation to be added]

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## Acknowledgments

- HTSlib developers for robust genomic file format handling
- OpenBLAS team for high-performance linear algebra routines
