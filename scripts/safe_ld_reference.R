#!/usr/bin/env Rscript
#
# safe_ld_reference.R - runnable reference implementation of SAFE-LD.
#
# This is a faithful port of the original `safe.ld.vcf()` R function, kept as the
# numerical reference that the C++ `safeld` tool is meant to reproduce. The
# algorithm is unchanged; what has been fixed is everything that made the
# original impossible to run somewhere else:
#
#   * hard-coded cluster paths for bcftools/plink2 are now arguments
#   * the MatrixEQTL / Matrix.utils / GEM dependencies are gone. MatrixEQTL was
#     only ever used as a chunked file reader, and the original already held the
#     whole dosage table in memory via fread() before handing it to MatrixEQTL,
#     so reading the matrix directly costs nothing and removes the dependency
#   * `seed` makes a run reproducible (the original could not be re-run)
#   * `W` can be loaded from a C++ run, which is what makes the two
#     implementations comparable at all (see "COMPARING" below)
#   * the intermediate .temp/.dose.gz round-trip through disk is dropped
#   * the hard-coded `##contig=<ID=29>` header line is replaced by the real
#     contigs read from the VCF
#   * `system()` calls are checked, and the temp dir is no longer removed with a
#     blind `rm -r`
#
# The arithmetic is deliberately identical to the original, including the
# details that look like quirks:
#   * scale() standardises each variant with the SAMPLE sd (n-1 denominator)
#   * the per-variant rescale to [0, 2] is min-max, applied across traits
#   * duplicate handling drops EVERY copy of a repeated ID, not just the extras
#
# COMPARING AGAINST THE C++ TOOL
# -----------------------------------------------------------------------------
# Both implementations draw the trait matrix W from rnorm(), so their outputs can
# never match cell-by-cell on separate runs - the randomness differs, not the
# method. To compare them properly, reuse the C++ tool's own W:
#
#   ./safeld preprocess -vcf in.vcf.gz -out prep -ntraits 100
#   ./safeld simulate   -prep prep -out results
#   ./safeld merge      -in results -out cpp.vcf.gz
#
#   Rscript safe_ld_reference.R --vcf in.vcf.gz --out r.vcf \
#           --ntraits 100 --w-from-prep prep
#
# `--w-from-prep` reads prep/traits/W_tile_*.bin, so both runs then use the same
# W and the dosages become directly comparable.
#
# Usage:
#   Rscript safe_ld_reference.R --vcf INPUT.vcf.gz --out OUTPUT.vcf [options]
#
#   --vcf FILE            input VCF (required)
#   --out FILE            output VCF (required)
#   --ntraits INT         number of synthetic traits (default 5000)
#   --maf FLOAT           MAF filter on INFO/AF (default 0.01)
#   --samples LIST        comma-separated sample IDs, or a file with one per line
#   --seed INT            RNG seed for the trait matrix
#   --slice-size INT      variants per slice (default 2000)
#   --w-from-prep DIR     load W from a C++ preprocess directory instead of rnorm
#   --bcftools PATH       bcftools binary (default: bcftools on PATH)
#   --impute-missing      mean-impute missing DS instead of erroring out
#   --dedup MODE          id (default, original R) | locus | both
#   --keep-tmp            keep the intermediate bcftools query output

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Extraction: the original vcf2matrixQtl(), minus the disk round-trip.
# ---------------------------------------------------------------------------

read_dosage_table <- function(vcf_file, bcftools = "bcftools", maf_filter = 0.01,
                              sample_list = NULL, dedup = "id",
                              impute_missing = FALSE, tmp_file = NULL,
                              keep_tmp = FALSE) {
  if (is.null(tmp_file)) {
    tmp_file <- tempfile(pattern = "safeld_query_", fileext = ".tsv")
  }
  if (!keep_tmp) {
    on.exit(unlink(tmp_file), add = TRUE)
  }

  # Same query as the original: ID, locus, alleles, INFO/AF and one DS per sample.
  fmt <- "%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF[\\t%DS]\\n"
  cmd <- sprintf("%s query -H -f '%s' %s > %s",
                 shQuote(bcftools), fmt, shQuote(vcf_file), shQuote(tmp_file))
  status <- system(cmd)
  if (status != 0) {
    stop("bcftools query failed (exit ", status, "). Command was:\n  ", cmd)
  }

  dosage <- fread(tmp_file, na.strings = c(".", "NA", ""), header = TRUE)
  if (nrow(dosage) == 0) {
    stop("bcftools query returned no records from ", vcf_file)
  }

  # Header comes back as "# [1]ID", "[2]CHROM", ..., "[7]SAMPLE:DS".
  header <- sub("^.*\\]", "", colnames(dosage))
  header <- sub(":DS$", "", header)
  header[1] <- "ID"
  colnames(dosage) <- header

  info_cols <- c("ID", "CHROM", "POS", "REF", "ALT", "AF")
  missing_cols <- setdiff(info_cols, colnames(dosage))
  if (length(missing_cols) > 0) {
    stop("unexpected bcftools output, missing column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  sample_cols <- setdiff(colnames(dosage), info_cols)
  if (length(sample_cols) == 0) {
    stop("no sample columns found; does the VCF carry a DS FORMAT field?")
  }

  # A multiallelic record makes bcftools emit "0.1,0.2" for AF, which is not
  # numeric. The original silently dropped these through the NA comparison
  # below; here it is at least reported.
  af <- suppressWarnings(as.numeric(dosage$AF))
  n_bad_af <- sum(is.na(af))
  if (n_bad_af > 0) {
    message("[R] dropping ", n_bad_af,
            " record(s) with missing or non-scalar INFO/AF (multiallelic?)")
  }
  dosage$AF <- af

  # MAF filter, exactly as the original: AF >= maf & AF <= 1 - maf.
  keep <- which(!is.na(dosage$AF) &
                dosage$AF >= maf_filter &
                dosage$AF <= (1 - maf_filter))
  message("[R] variants after MAF filter: ", length(keep), "/", nrow(dosage))
  dosage <- dosage[keep, ]
  if (nrow(dosage) == 0) {
    stop("no variants passed the MAF filter")
  }

  # Duplicate handling. The original keyed on ID alone and removed every copy.
  # That is safe on the imputed VCFs it was written for, where %ID is always
  # populated, but it silently discards the entire file when the ID column is
  # "." - which is why the C++ tool now keys on the locus. Mode is selectable so
  # the two can be lined up.
  dup_key <- switch(dedup,
    id    = dosage$ID,
    locus = paste(dosage$CHROM, dosage$POS, dosage$REF, dosage$ALT, sep = ":"),
    both  = NULL,
    stop("--dedup must be one of: id, locus, both"))

  if (identical(dedup, "both")) {
    locus <- paste(dosage$CHROM, dosage$POS, dosage$REF, dosage$ALT, sep = ":")
    dup <- duplicated(locus) | duplicated(locus, fromLast = TRUE)
    has_id <- !is.na(dosage$ID) & dosage$ID != "."
    dup_id <- has_id & (duplicated(dosage$ID) | duplicated(dosage$ID, fromLast = TRUE))
    drop <- dup | dup_id
  } else {
    if (identical(dedup, "id")) {
      n_placeholder <- sum(is.na(dosage$ID) | dosage$ID == ".")
      if (n_placeholder > 0) {
        warning("--dedup id: ", n_placeholder, " variant(s) have no ID ('.'). ",
                "They all share one key and will ALL be dropped as duplicates. ",
                "Use --dedup locus for such a VCF.", call. = FALSE)
      }
    }
    drop <- duplicated(dup_key) | duplicated(dup_key, fromLast = TRUE)
  }

  if (any(drop)) {
    message("[R] dropping ", sum(drop), " variant(s) as duplicates (mode: ", dedup, ")")
    dosage <- dosage[!drop, ]
  }
  if (nrow(dosage) == 0) {
    stop("every variant was removed by deduplication")
  }

  # Sample subsetting, as in the original: INFO/AF above still describes the FULL
  # cohort, so the MAF filter is not recomputed for the subset.
  if (!is.null(sample_list)) {
    found <- intersect(sample_list, sample_cols)
    if (length(found) == 0) {
      stop("none of the requested samples are present in the VCF")
    }
    if (length(found) < length(sample_list)) {
      message("[R] ", length(found), "/", length(sample_list),
              " requested samples found in the VCF")
    }
    warning("sample subsetting keeps the whole-cohort INFO/AF for MAF filtering, ",
            "matching the original R code", call. = FALSE)
    sample_cols <- found
  }

  G <- as.matrix(dosage[, ..sample_cols])
  storage.mode(G) <- "double"

  n_missing <- sum(is.na(G))
  if (n_missing > 0) {
    if (!impute_missing) {
      stop(n_missing, " missing DS value(s) found. The original R code has no ",
           "missingness handling and would produce NA dosages here. Re-run with ",
           "--impute-missing to mean-impute (what the C++ tool does).")
    }
    # Mean of the observed calls, per variant - matches the C++ imputation.
    row_means <- rowMeans(G, na.rm = TRUE)
    idx <- which(is.na(G), arr.ind = TRUE)
    G[idx] <- row_means[idx[, "row"]]
    message("[R] mean-imputed ", n_missing, " missing DS value(s)")
  }

  list(info = dosage[, ..info_cols], G = G, samples = sample_cols)
}

# ---------------------------------------------------------------------------
# Trait matrix
# ---------------------------------------------------------------------------

# Reproduces `matrice` from the original: ntraits x nsamples of standard normals.
make_traits <- function(ntraits, nsamples, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  matrix(rnorm(ntraits * nsamples), nrow = ntraits, ncol = nsamples)
}

# Reads the trait matrix written by the C++ preprocess stage, so both
# implementations can be run against identical W.
load_w_from_prep <- function(prep_dir) {
  meta_file <- file.path(prep_dir, "traits", "metadata.txt")
  if (!file.exists(meta_file)) {
    stop("no trait metadata at ", meta_file)
  }
  lines <- readLines(meta_file)
  kv <- regmatches(lines, regexpr("=", lines), invert = TRUE)
  keys <- vapply(kv, `[`, character(1), 1)
  vals <- vapply(kv, `[`, character(1), 2)
  get <- function(k) {
    i <- match(k, keys)
    if (is.na(i)) stop("missing key '", k, "' in ", meta_file)
    vals[i]
  }

  n_traits  <- as.integer(get("n_traits"))
  n_samples <- as.integer(get("n_samples"))
  counts    <- as.integer(strsplit(get("tile_trait_counts"), ",")[[1]])

  W <- matrix(0, nrow = n_traits, ncol = n_samples)
  row0 <- 0L
  for (i in seq_along(counts)) {
    tile_file <- file.path(prep_dir, "traits", sprintf("W_tile_%d.bin", i - 1L))
    con <- file(tile_file, "rb")
    v <- readBin(con, "double", n = counts[i] * n_samples, size = 8)
    close(con)
    if (length(v) != counts[i] * n_samples) {
      stop("short read on ", tile_file)
    }
    # Tiles are row-major traits x samples.
    W[(row0 + 1L):(row0 + counts[i]), ] <- matrix(v, nrow = counts[i],
                                                  ncol = n_samples, byrow = TRUE)
    row0 <- row0 + counts[i]
  }
  message("[R] loaded W from ", prep_dir, ": ", n_traits, " traits x ",
          n_samples, " samples")
  W
}

# ---------------------------------------------------------------------------
# Core arithmetic
# ---------------------------------------------------------------------------

# Mirrors `apply(slice_mat, MARGIN = 1, scale)`: standardises each variant (row)
# across samples and returns the transpose, i.e. samples x variants.
# scale() divides by the SAMPLE sd (n-1 denominator); kept as-is.
scale_variants <- function(m) {
  centered <- m - rowMeans(m)
  sdev <- sqrt(rowSums(centered^2) / (ncol(m) - 1))
  t(centered / sdev)
}

# Mirrors `dosager`: shift to zero, divide by the max, scale to [0, 2].
# Equivalent to 2 * (x - min) / (max - min).
dosager <- function(x) {
  x <- x - min(x)
  x <- x / max(x)
  x * 2
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

safe_ld_vcf <- function(vcf_file, out_file, ntraits = 5000, maf_filter = 0.01,
                        sample_list = NULL, seed = NULL, slice_size = 2000,
                        bcftools = "bcftools", w_from_prep = NULL,
                        impute_missing = FALSE, dedup = "id", keep_tmp = FALSE) {
  t0 <- Sys.time()
  message("[R] SAFE-LD reference implementation, started ", format(t0))

  dat <- read_dosage_table(vcf_file, bcftools = bcftools, maf_filter = maf_filter,
                           sample_list = sample_list, dedup = dedup,
                           impute_missing = impute_missing, keep_tmp = keep_tmp)
  info <- dat$info
  G <- dat$G
  n_variants <- nrow(G)
  n_samples <- ncol(G)
  message("[R] genotype matrix: ", n_variants, " variants x ", n_samples, " samples")

  if (n_samples < 2) {
    stop("at least 2 samples are required to standardise genotypes")
  }

  if (!is.null(w_from_prep)) {
    W <- load_w_from_prep(w_from_prep)
    if (ncol(W) != n_samples) {
      stop("W has ", ncol(W), " samples but the genotype matrix has ", n_samples,
           ". The C++ run and this run must use the same sample set.")
    }
    ntraits <- nrow(W)
  } else {
    W <- make_traits(ntraits, n_samples, seed = seed)
    message("[R] generated W: ", ntraits, " traits x ", n_samples, " samples",
            if (is.null(seed)) " (no seed: not reproducible)" else
              paste0(" (seed ", seed, ")"))
  }

  if (ntraits < 2) {
    warning("with fewer than 2 traits the per-variant min-max rescale has no ",
            "spread to work with and every dosage collapses to the same value",
            call. = FALSE)
  }

  # Zero-variance variants make scale() produce NaN, which the original then
  # wrote straight into the VCF. Report them rather than letting them pass silently.
  variant_sd <- sqrt(rowSums((G - rowMeans(G))^2) / (n_samples - 1))
  n_constant <- sum(!(variant_sd > 0))
  if (n_constant > 0) {
    warning(n_constant, " variant(s) have zero variance across the selected ",
            "samples and will produce NaN dosages (the C++ tool drops them)",
            call. = FALSE)
  }

  out <- matrix(NA_real_, nrow = n_variants, ncol = ntraits,
                dimnames = list(NULL, paste0("T", seq_len(ntraits))))

  n_slices <- max(1L, ceiling(n_variants / slice_size))
  for (i in seq_len(n_slices)) {
    lo <- (i - 1L) * slice_size + 1L
    hi <- min(i * slice_size, n_variants)

    # samples x variants, standardised per variant
    slice_mat <- scale_variants(G[lo:hi, , drop = FALSE])

    # (traits x samples) %*% (samples x variants) -> traits x variants
    betas <- W %*% slice_mat
    betas <- t(betas)              # variants x traits
    betas <- betas / n_samples

    # rescale each variant across traits into [0, 2]
    betas <- t(apply(betas, MARGIN = 1, dosager))

    out[lo:hi, ] <- betas
    message("[R] slice ", i, "/", n_slices, " (variants ", lo, "-", hi, ")")
  }

  write_vcf(out_file, info, out, vcf_file, bcftools)

  t1 <- Sys.time()
  message("[R] SAFE-LD complete, wrote ", out_file, " in ",
          format(round(difftime(t1, t0), 2)))
  invisible(out)
}

# The original wrote a hard-coded `##contig=<ID=29>`; take the real contigs from
# the input header instead so the output can be indexed.
write_vcf <- function(out_file, info, dosages, vcf_file, bcftools) {
  contigs <- character(0)
  hdr <- suppressWarnings(system(sprintf("%s view -h %s", shQuote(bcftools),
                                         shQuote(vcf_file)), intern = TRUE))
  if (!is.null(attr(hdr, "status")) && attr(hdr, "status") != 0) {
    warning("could not read the VCF header for contig lines", call. = FALSE)
  } else {
    contigs <- grep("^##contig=", hdr, value = TRUE)
  }

  header <- c(
    "##fileformat=VCFv4.1",
    "##source=safe_ld_reference.R",
    contigs,
    "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Dosage\">"
  )
  writeLines(header, out_file)

  vcf <- data.table(
    `#CHROM` = info$CHROM,
    POS      = info$POS,
    ID       = info$ID,
    REF      = info$REF,
    ALT      = info$ALT,
    QUAL     = ".",
    FILTER   = "PASS",
    INFO     = ".",
    FORMAT   = "DS"
  )
  vcf <- cbind(vcf, as.data.table(round(dosages, 4)))

  fwrite(vcf, file = out_file, sep = "\t", quote = FALSE,
         row.names = FALSE, col.names = TRUE, append = TRUE)
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  opts <- list(vcf = NULL, out = NULL, ntraits = 5000, maf = 0.01,
               samples = NULL, seed = NULL, slice_size = 2000,
               bcftools = "bcftools", w_from_prep = NULL,
               impute_missing = FALSE, dedup = "id", keep_tmp = FALSE)

  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    take <- function() {
      if (i + 1L > length(args)) stop(a, " needs a value")
      args[i + 1L]
    }
    switch(a,
      "--vcf"            = { opts$vcf <- take(); i <- i + 1L },
      "--out"            = { opts$out <- take(); i <- i + 1L },
      "--ntraits"        = { opts$ntraits <- as.integer(take()); i <- i + 1L },
      "--maf"            = { opts$maf <- as.numeric(take()); i <- i + 1L },
      "--samples"        = { opts$samples <- take(); i <- i + 1L },
      "--seed"           = { opts$seed <- as.integer(take()); i <- i + 1L },
      "--slice-size"     = { opts$slice_size <- as.integer(take()); i <- i + 1L },
      "--bcftools"       = { opts$bcftools <- take(); i <- i + 1L },
      "--w-from-prep"    = { opts$w_from_prep <- take(); i <- i + 1L },
      "--dedup"          = { opts$dedup <- take(); i <- i + 1L },
      "--impute-missing" = { opts$impute_missing <- TRUE },
      "--keep-tmp"       = { opts$keep_tmp <- TRUE },
      "--help"           = { opts$help <- TRUE },
      "-h"               = { opts$help <- TRUE },
      stop("unknown argument: ", a))
    i <- i + 1L
  }
  opts
}

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    cat(paste(readLines(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), n = 60),
        collapse = "\n"), "\n")
    quit(status = if (length(args) == 0L) 1L else 0L)
  }

  opts <- parse_args(args)
  if (is.null(opts$vcf) || is.null(opts$out)) {
    stop("--vcf and --out are required")
  }

  sample_list <- NULL
  if (!is.null(opts$samples)) {
    sample_list <- if (file.exists(opts$samples)) {
      trimws(readLines(opts$samples))
    } else {
      trimws(strsplit(opts$samples, ",")[[1]])
    }
    sample_list <- sample_list[nzchar(sample_list)]
  }

  safe_ld_vcf(vcf_file = opts$vcf, out_file = opts$out, ntraits = opts$ntraits,
              maf_filter = opts$maf, sample_list = sample_list, seed = opts$seed,
              slice_size = opts$slice_size, bcftools = opts$bcftools,
              w_from_prep = opts$w_from_prep, impute_missing = opts$impute_missing,
              dedup = opts$dedup, keep_tmp = opts$keep_tmp)
}
