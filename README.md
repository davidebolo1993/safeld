# SAFELD

A high-performance C++ implementation of the SAFE-LD (Shrinkage and Anonymization Framework for LD Estimation) method for generating synthetic genomic data from VCF files.

## Overview

SAFELD processes VCF files to generate synthetic traits while preserving the linkage disequilibrium structure of the original data. This tool is useful for:

- Generating synthetic genomic datasets for testing and validation
- Privacy-preserving genomic data sharing
- Method development and benchmarking in genomics research

## Features

- **High Performance**: Optimized C++ implementation with OpenBLAS integration and BLAS GEMM operations
- **Memory Efficient**: Streaming, chunked data structures that keep the working set bounded
- **Parallel Processing**: Multi-threaded BLAS, with the thread count exposed via `-workers`
- **Flexible Input**: Supports compressed and uncompressed VCF files
- **HTSlib Integration**: Robust VCF parsing using industry-standard library
- **Docker Support**: Containerized deployment for reproducibility
- **Scalable Workflow**: Three-stage pipeline (`preprocess`, `simulate`, `merge`) for large cohorts and high trait counts
- **Tile Streaming**: Trait tiles are generated and consumed on demand to keep RAM bounded
- **Configurable Batching**: Simulation variant batch size is tunable via CLI

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
- `safeld`

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
docker run --rm -v $(pwd):/data safeld preprocess -vcf /data/input.vcf.gz -out /data/prep
```

## Usage

### Three-Stage Workflow

SAFELD uses a **three-stage workflow**:

Global option:
- `-verbose` Enable detailed debug logging

#### Stage 1: Preprocessing

Generates trait matrix and partitions variants into chunks. Accepts a VCF, a
plink2 `.pgen/.pvar/.psam` triple, or a plink1 `.bed/.bim/.fam` triple — exactly
one of `-vcf`, `-pfile`, `-bfile`. The plink formats are read natively, with no
VCF conversion, when built with `-DSAFELD_PGEN=ON` (see below).
Input VCF must be coordinate-sorted and biallelic (split multiallelic records
first with `bcftools norm -m -any`; non-biallelic records are skipped and counted).

```bash
./safeld preprocess \
  -vcf input.vcf.gz \
  -out preprocessed_data \
  -ntraits 5000 \
  -maf 0.01
```

**Preprocessing Options:**

```txt
./safeld preprocess [OPTIONS]

Options:
  -vcf FILE            Input VCF file
  -pfile PREFIX        Input plink2 .pgen/.pvar/.psam
  -bfile PREFIX        Input plink1 .bed/.bim/.fam
  -out DIR             Output directory for preprocessed data
  -samples LIST        Comma-separated sample IDs, or a file with one per line
  -extract FILE        Keep only these variant IDs, one per line
  -maf FLOAT           MAF filter (default: 0.01)
  -max-missing FLOAT   Max fraction of missing calls per variant (default: 0.1)
  -dosage-field FIELD  auto|DS|GT: which FORMAT field to read (default: auto)
                       auto measures both fields up front and picks; DS fills
                       absent dosages from GT; GT ignores DS entirely
  -ntraits INT         Number of traits (default: 10)
  -chunk-size INT      Variants per chunk (default: 10000)
  -traits-per-tile INT Traits per tile (default: auto, ~1GB tiles)
  -h, --help           Show this help message
```

> **Disk space:** preprocessing writes a temporary deduplication spool into the
> output directory holding one copy of the genotype matrix
> (`n_variants x n_samples x 8` bytes). It is removed when the stage finishes.

**Output structure:**
```
preprocessed_data/
├── header_contigs.txt    # Serialized ##contig header lines
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
  -variant-batch-size INT Variants per simulation batch (default: 4000)
  -compress          Compress output VCF chunks
  -start-chunk INT   First chunk to process (default: all)
  -end-chunk INT     Last chunk to process (default: all)
  -h, --help         Show this help message
```


#### Stage 3: Merge

Combines all chunk VCF files into a single output file.
When output is compressed, a `.tbi` index is created by default.
If merged chunks are unexpectedly unsorted, merge falls back to `bcftools sort`.

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

## Algorithm

SAFELD separates preprocessing from simulation for scalability:

**Preprocessing Stage:**
1. Parse VCF once, skip non-biallelic records, and apply MAF and missingness filtering
2. Deduplicate on locus (`CHROM:POS:REF:ALT`) and, where present, on variant ID
3. Generate trait matrix W (T × S) with standard normal random values
4. Standardize genotype dosages and partition into chunks (B variants each)
5. Serialize traits (tiled if large) and genotype chunks to disk

**Missing genotypes:**

Missingness is resolved during preprocessing, before standardization:

- A `DS` value that is absent, negative or `NaN` counts as missing. Without `DS`,
  dosages are derived from `GT` on the diploid 0–2 scale, normalized by each
  sample's own ploidy so a hemizygous ALT call (chrX/chrY in a male,
  mitochondria) scores 2.0 rather than being confused with a heterozygote.
- **`preprocess` scans the head of the input before it starts** and reports what
  it found: sample count, non-biallelic records, and the fraction of calls
  carrying `GT` and `DS`. In `auto` mode it then picks the dosage source from
  those measurements and says which it chose and why, so a run explains its own
  input without a separate diagnostic step.
- **`-dosage-field` matters for VCFs that carry both `GT` and `DS`.** Some
  exports (plink2 in particular) write the `DS` subfield for only a fraction of
  samples while `GT` stays complete — `0|1:0.97` sitting next to a bare `0|0`.
  Such a sample is *not* missing: its genotype is known from `GT`. The default
  `auto` therefore fills each absent dosage from that sample's own hard call
  rather than imputing it, and reports how many calls it filled. Treating those
  gaps as missing and mean-imputing them attenuates every pairwise r² in
  proportion to the `DS` presence rate — on a file with 30% `DS` coverage a true
  r² of 0.80 reads as 0.09 — which appears as distinct lower bands against the
  original LD. `-dosage-field GT` builds the matrix from hard calls throughout;
  `-dosage-field DS` keeps whatever `DS` exists and fills the rest from `GT`.
- A genotype with **any** missing allele (`./1`) counts as missing outright; it
  is not silently scored as a reference call.
- Variants whose missing fraction exceeds `-max-missing` are dropped and counted.
- Surviving missing entries are replaced with the mean of that variant's observed
  dosages, which leaves the variant mean unchanged. Because imputed entries sit
  exactly at the mean, σ is shrunk by `sqrt(n_observed / n)`; this is the usual
  convention (plink does the same).
- Allele frequency and missingness are measured on observed calls only. `INFO/AF`
  is used as a fast pre-filter only when the whole cohort is in use — under
  `-samples` it describes samples that are not in the matrix, so AF is recomputed
  from the selected samples instead.
- Variants with no variance across the selected samples are dropped, since
  per-variant scaling would otherwise emit them as a constant dosage for every trait.

**Simulation Stage:**
1. Load trait metadata once
2. For each chunk:
   - Load standardized genotypes G (B × S)
   - Process variants in configurable batches (`-variant-batch-size`)
   - Stream trait tiles and compute Y = G × W^T using optimized BLAS GEMM
   - Scale synthetic dosages to [0, 2] range
   - Write chunk VCF to disk
3. Release chunk memory before processing next

**Merge Stage:**
- Concatenate all chunk VCF files maintaining chromosome order
- Stream processing for memory efficiency


### Native .pgen / .bed support

Reading plink formats directly avoids the VCF round-trip that caused real
trouble: exporting a hard-call pgen with `--export vcf vcf-dosage=DS` writes the
`DS` subfield for only the entries that happen to carry a dosage track, and the
resulting VCF cannot be read correctly without knowing that. A `.pgen` records
dosage presence explicitly per sample, so the ambiguity does not arise.

```bash
git clone --depth 1 https://github.com/chrchang/plink-ng
cmake -DSAFELD_PGEN=ON -DPLINK_NG_DIR=/path/to/plink-ng ..
```

A full clone is needed: plink-ng's vendored `simde/` is required to compile.
pgenlib is LGPL-3.0 and is built as a shared library; safeld itself stays MIT.
Without this option `-pfile`/`-bfile` report a clear error and `-vcf` is
unaffected.

### Subsetting

`-extract FILE` keeps only the listed variant IDs, matching the ID column of the
VCF or `.pvar`/`.bim` — the same semantics as plink's `--extract`:

```txt
1:113989901:A:G
1:113990655:A:G
```

`-samples` takes either a comma-separated list or a file with one ID per line.
Both work for every input format. If entries in the extract list match nothing,
preprocessing says so rather than quietly keeping fewer variants than expected.

## Tests

```bash
tests/run_tests.sh build/safeld
```

The central check writes the same genotypes as both a VCF and a plink1 `.bed`,
then requires the two readers to produce byte-identical standardized matrices —
so a reader bug shows up as disagreement rather than as plausible output. The
rest cover the input pathologies that caused real failures: ID-less records,
duplicate loci and IDs, multiallelic sites, all-missing and monomorphic
variants, and sparse `DS`. Tests needing `.bed` skip themselves when built
without `SAFELD_PGEN`.

### Dependencies

- **HTSlib**: VCF/BCF file format handling
- **OpenBLAS**: Optimized linear algebra operations (cblas_ddot, cblas_dgemm)
- **OpenMP**: Parallel processing support
- **BGZF**: Block gzip compression for output
- **bcftools**: Optional fallback sorting if merged chunks are unexpectedly out of order

## File Formats

### Input VCF Requirements

- Must contain `DS` (dosage) format field or `GT` (genotype) field
- Should include `AF` (allele frequency) in INFO field (calculated if missing, and
  always recomputed from the selected samples when `-samples` is used)
- Supports both compressed (.vcf.gz) and uncompressed (.vcf) files
- Must be coordinate-sorted (enforced during preprocessing)
- Must be biallelic; multiallelic records are skipped and reported. Split them
  with `bcftools norm -m -any` to keep them

### Output VCF Structure

```txt
##fileformat=VCFv4.1
##source=safeld
##contig=<ID=22>
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

If you use SAFELD in your research, please cite our [preprint](https://doi.org/10.1101/2025.09.29.679154).

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## Acknowledgments

- HTSlib developers for robust genomic file format handling
- OpenBLAS team for high-performance linear algebra routines
