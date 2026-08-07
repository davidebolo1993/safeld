#!/usr/bin/env bash
#
# ld_check.sh - run safeld two ways from one .pgen and report how well the
# simulated LD matches the original.
#
#   route A   plink2 exports a pre-extracted VCF, safeld reads that
#   route B   safeld reads the .pgen directly and does the extraction itself
#
# Both should land on the same answer. Each route gets its own output directory
# so nothing overwrites anything.
#
# Usage:
#   scripts/ld_check.sh -p PGEN_PREFIX -e EXTRACT_FILE -o OUTDIR [options]
#
#   -p PREFIX     .pgen/.pvar/.psam prefix (required)
#   -e FILE       variant IDs to keep, one per line (required)
#   -o DIR        output directory (required)
#   -n INT        traits (default 10000)
#   -m FLOAT      MAF threshold (default 0.01)
#   -w INT        LD window in kb (default 1000)
#   -s PATH       safeld binary (default: build/safeld next to this script)
#   -k PATH       plink2 binary (default: plink2 from PATH)
#
# Example:
#   scripts/ld_check.sh -p /path/GRCh38_ukb_processed_autosomes \
#                       -e ../snp_chr1.txt -o ldcheck_chr1

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

PGEN=""; EXTRACT=""; OUT=""
NTRAITS=10000; MAF=0.01; WINDOW=1000
SAFELD="$HERE/../build/safeld"
PLINK2="plink2"

while getopts "p:e:o:n:m:w:s:k:h" opt; do
  case "$opt" in
    p) PGEN="$OPTARG" ;;
    e) EXTRACT="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    n) NTRAITS="$OPTARG" ;;
    m) MAF="$OPTARG" ;;
    w) WINDOW="$OPTARG" ;;
    s) SAFELD="$OPTARG" ;;
    k) PLINK2="$OPTARG" ;;
    h) sed -n '2,26p' "$0"; exit 0 ;;
    *) exit 2 ;;
  esac
done

[ -n "$PGEN" ] && [ -n "$EXTRACT" ] && [ -n "$OUT" ] || {
  echo "need -p, -e and -o; see -h" >&2; exit 2; }

PGEN="${PGEN%.pgen}"

mkdir -p "$OUT"/{input,from_vcf,from_pgen,original}
OUT="$(cd "$OUT" && pwd)"
EXTRACT="$(cd "$(dirname "$EXTRACT")" && pwd)/$(basename "$EXTRACT")"

echo "safeld  : $SAFELD"
echo "plink2  : $PLINK2"
echo "pgen    : $PGEN"
echo "extract : $EXTRACT"
echo "traits  : $NTRAITS   maf: $MAF   window: ${WINDOW}kb"
echo "output  : $OUT"
echo

# ---------------------------------------------------------------------------
# 1. Pre-extracted VCF, the way it was done by hand. DS-force writes a dosage
#    for every sample, substituting the hard call where none is stored.
# ---------------------------------------------------------------------------
echo "[1/6] exporting pre-extracted VCF"
"$PLINK2" --pfile "$PGEN" --extract "$EXTRACT" --maf "$MAF" \
          --export vcf vcf-dosage=DS-force \
          --out "$OUT/input/extracted" > "$OUT/input/export.log" 2>&1

# ---------------------------------------------------------------------------
# 2. Route A: safeld reads that VCF.
# ---------------------------------------------------------------------------
echo "[2/6] safeld from the VCF"
"$SAFELD" preprocess -vcf "$OUT/input/extracted.vcf" -maf "$MAF" \
          -ntraits "$NTRAITS" -out "$OUT/from_vcf/prep"      > "$OUT/from_vcf/preprocess.log" 2>&1
"$SAFELD" simulate -prep "$OUT/from_vcf/prep" -compress \
          -out "$OUT/from_vcf/sim"                            > "$OUT/from_vcf/simulate.log" 2>&1
"$SAFELD" merge -in "$OUT/from_vcf/sim" \
          -out "$OUT/from_vcf/safeld.vcf.gz"                  > "$OUT/from_vcf/merge.log" 2>&1

# ---------------------------------------------------------------------------
# 3. Route B: safeld reads the whole .pgen and extracts internally.
# ---------------------------------------------------------------------------
echo "[3/6] safeld from the pgen"
"$SAFELD" preprocess -pfile "$PGEN" -extract "$EXTRACT" -maf "$MAF" \
          -ntraits "$NTRAITS" -out "$OUT/from_pgen/prep"     > "$OUT/from_pgen/preprocess.log" 2>&1
"$SAFELD" simulate -prep "$OUT/from_pgen/prep" -compress \
          -out "$OUT/from_pgen/sim"                           > "$OUT/from_pgen/simulate.log" 2>&1
"$SAFELD" merge -in "$OUT/from_pgen/sim" \
          -out "$OUT/from_pgen/safeld.vcf.gz"                 > "$OUT/from_pgen/merge.log" 2>&1

# ---------------------------------------------------------------------------
# 4. LD of the real genotypes.
# ---------------------------------------------------------------------------
echo "[4/6] LD of the original genotypes"
"$PLINK2" --vcf "$OUT/input/extracted.vcf" dosage=DS --make-pgen \
          --out "$OUT/original/pgen"      > "$OUT/original/import.log" 2>&1
"$PLINK2" --pfile "$OUT/original/pgen" --r2-unphased \
          --ld-window-kb "$WINDOW" --ld-window-r2 0 \
          --out "$OUT/original/ld"        > "$OUT/original/ld.log" 2>&1

# ---------------------------------------------------------------------------
# 5. LD of each simulated set.
# ---------------------------------------------------------------------------
echo "[5/6] LD of the simulated sets"
for route in from_vcf from_pgen; do
  "$PLINK2" --vcf "$OUT/$route/safeld.vcf.gz" dosage=DS --make-pgen \
            --out "$OUT/$route/pgen"     > "$OUT/$route/import.log" 2>&1
  "$PLINK2" --pfile "$OUT/$route/pgen" --r2-unphased \
            --ld-window-kb "$WINDOW" --ld-window-r2 0 \
            --out "$OUT/$route/ld"       > "$OUT/$route/ld.log" 2>&1
done

# ---------------------------------------------------------------------------
# 6. Correlate. ld_cor.awk joins pairs on position rather than assuming the two
#    files list the same pairs in the same order.
# ---------------------------------------------------------------------------
echo "[6/6] correlating"
{
  echo "safeld LD check"
  echo "  pgen    : $PGEN"
  echo "  extract : $EXTRACT"
  echo "  traits  : $NTRAITS   maf: $MAF   window: ${WINDOW}kb"
  echo "  date    : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  for route in from_vcf from_pgen; do
    echo
    echo "=== $route ==="
    grep -E 'kept|excluded' "$OUT/$route/preprocess.log" || true
    awk -f "$HERE/ld_cor.awk" "$OUT/original/ld.vcor" "$OUT/$route/ld.vcor"
  done
} | tee "$OUT/summary.txt"

echo
echo "summary written to $OUT/summary.txt"
