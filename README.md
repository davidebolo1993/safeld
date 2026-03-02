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
- **Scalable Workflow**: Three-stage pipeline (`preprocess`, `simulate`, `merge`) for large cohorts and high trait counts
- **Tile Streaming**: Trait tiles are generated and consumed on demand to keep RAM bounded
- **VCF Compliance**: Chunk outputs include `##contig` header records for downstream tools
- **Sorted Merge + Indexing**: Merge validates sort order and builds tabix indexes in-house (HTSlib)

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

This builds one executable:
- `safeld` - Three-stage workflow (memory-efficient and scalable)

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

# Run with Docker (example: stage 1 preprocess)
docker run --rm -v $(pwd):/data safeld preprocess -vcf /data/input.vcf.gz -out /data/prep -ntraits 10000
```

## Usage

### Three-Stage Workflow

SAFELD uses a **three-stage workflow**:

#### Stage 1: Preprocessing

Generates trait matrix and partitions variants into chunks.
Input VCF must be coordinate-sorted.

```bash
./safeld preprocess \
  -vcf input.vcf.gz \
  -out preprocessed_data \
  -ntraits 10000 \
  -chunk-size 5000 \
  -maf 0.01
```

**Preprocessing Options:**

```txt
./safeld preprocess [OPTIONS]

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
./safeld simulate \
  -prep preprocessed_data \
  -out results \
  -workers 32 \
  -compress

# Process specific chunk range (useful for cluster parallelization)
./safeld simulate \
  -prep preprocessed_data \
  -out results \
  -start-chunk 0 \
  -end-chunk 9 \
  -workers 32 \
  -compress
```

**Simulation Options:**

```txt
./safeld simulate [OPTIONS]

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
./safeld merge \
  -in results \
  -out final_output.vcf.gz
```

**Merge Options:**

```txt
./safeld merge [OPTIONS]

Options:
  -in DIR            Directory with chunk VCF files (required)
  -out FILE          Output merged VCF file (required)
  -no-compress       Don't compress output (default: compressed)
  -no-index          Don't create tabix index for compressed output
  -no-sort           Skip sortedness enforcement during merge
  -h, --help         Show this help message
```

### Complete Workflow Example

```bash
# 1. Preprocess (one time per dataset)
./safeld preprocess \
  -vcf large_dataset.vcf.gz \
  -out preprocessed_chr22 \
  -ntraits 10000 \
  -chunk-size 5000 \
  -maf 0.01

# 2. Simulate (can be parallelized)
./safeld simulate \
  -prep preprocessed_chr22 \
  -out results_chr22 \
  -workers 32 \
  -compress

# 3. Merge
./safeld merge \
  -in results_chr22 \
  -out chr22_final.vcf.gz
```

## Algorithm

SAFELD separates preprocessing from simulation for scalability:

**Preprocessing Stage:**
1. Parse VCF once and apply MAF filtering
2. Generate trait matrix W (T × S) with standard normal random values
3. Standardize genotype dosages and partition into chunks (B variants each)
4. Serialize traits (tiled if large) and genotype chunks to disk

**Simulation Stage:**
1. Load trait metadata once
2. For each chunk:
   - Load standardized genotypes G (B × S)
   - Stream trait tiles and compute Y = G × W^T using optimized BLAS GEMM
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
- **bcftools**: Optional fallback sorting if merged chunks are unexpectedly out of order

## File Formats

### Input VCF Requirements

- Must contain `DS` (dosage) format field or `GT` (genotype) field
- Should include `AF` (allele frequency) in INFO field (calculated if missing)
- Supports both compressed (.vcf.gz) and uncompressed (.vcf) files
- Must be coordinate-sorted (enforced during preprocessing)

### Output VCF Structure

```txt
##fileformat=VCFv4.1
##source=safeld
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT  T1    T2    ...
chr1    1000    rs123   A    G    .     PASS    .     DS      1.23  0.45  ...
```

### Binary Formats

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
