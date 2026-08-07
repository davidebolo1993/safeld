#!/usr/bin/env bash
#
# run_tests.sh - end-to-end checks for safeld.
#
# The central test is a cross-check rather than a golden file: the same
# genotypes are written as a VCF and as a plink1 .bed, and the two readers must
# produce byte-identical standardized matrices. A golden file only tells you
# that today's output matches yesterday's; this tells you that two independent
# code paths agree on what the data means.
#
# Usage:
#   tests/run_tests.sh [path/to/safeld] [workdir]
#
# Tests needing .bed/.pgen are skipped automatically when the binary was built
# without SAFELD_PGEN.

set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SAFELD="${1:-$HERE/../build/safeld}"
WORK="${2:-$HERE/.work}"
REAL="${3:-$HERE/realdata}"      # optional; see get_real_testdata.sh

if [ ! -x "$SAFELD" ]; then
  echo "safeld binary not found at $SAFELD" >&2
  echo "usage: tests/run_tests.sh [path/to/safeld] [workdir]" >&2
  exit 2
fi

PASS=0; FAIL=0; SKIP=0
ok()   { printf "  \033[32mPASS\033[0m  %s\n" "$1"; PASS=$((PASS+1)); }
bad()  { printf "  \033[31mFAIL\033[0m  %s\n" "$1"; [ $# -gt 1 ] && printf "        %s\n" "$2"; FAIL=$((FAIL+1)); }
skip() { printf "  \033[33mSKIP\033[0m  %s\n" "$1"; SKIP=$((SKIP+1)); }

rm -rf "$WORK"; mkdir -p "$WORK"
DATA="$WORK/data"
python3 "$HERE/make_testdata.py" "$DATA" >/dev/null || { echo "fixture generation failed" >&2; exit 2; }

# Does this build have .pgen/.bed support? Decide on the message, not the exit
# code, so the check does not depend on how failures are reported.
probe=$("$SAFELD" preprocess -bfile "$DATA/matched" -out "$WORK/_probe" -ntraits 2 2>&1 || true)
if grep -q "no .pgen/.bed support" <<<"$probe"; then HAVE_PGEN=0; else HAVE_PGEN=1; fi
rm -rf "$WORK/_probe"

run() { "$SAFELD" "$@" 2>&1; }
kept() { grep -oE 'kept [0-9,]+' <<<"$1" | tail -1 | tr -d 'kept ,'; }

echo
echo "safeld tests  ($SAFELD)"
echo

# --------------------------------------------------------------- readers ----
echo "readers"
out=$(run preprocess -vcf "$DATA/matched.vcf" -out "$WORK/a_vcf" -ntraits 20)
[ "$(kept "$out")" = "60" ] && ok "VCF reader keeps all 60 clean variants" \
                            || bad "VCF reader keeps all 60 clean variants" "kept $(kept "$out")"

if [ "$HAVE_PGEN" = "1" ]; then
  out=$(run preprocess -bfile "$DATA/matched" -out "$WORK/a_bed" -ntraits 20)
  [ "$(kept "$out")" = "60" ] && ok "bed/bim/fam reader keeps all 60" \
                              || bad "bed/bim/fam reader keeps all 60" "kept $(kept "$out")"

  # The point of the whole fixture: two readers, one answer.
  if cmp -s "$WORK/a_vcf/chunks/chunk_0.bin" "$WORK/a_bed/chunks/chunk_0.bin"; then
    ok "VCF and bed readers agree byte for byte"
  else
    bad "VCF and bed readers agree byte for byte" "chunk_0.bin differs"
  fi
  if diff -q <(grep '^1' "$WORK/a_vcf/chunks/chunk_0.meta") \
             <(grep '^1' "$WORK/a_bed/chunks/chunk_0.meta") >/dev/null; then
    ok "VCF and bed readers agree on variant metadata"
  else
    bad "VCF and bed readers agree on variant metadata"
  fi
else
  skip "bed/bim/fam reader (built without SAFELD_PGEN)"
  skip "VCF vs bed byte equality (built without SAFELD_PGEN)"
fi

# -------------------------------------------------------------- subsetting --
echo
echo "subsetting"
out=$(run preprocess -vcf "$DATA/matched.vcf" -out "$WORK/b_vcf" -ntraits 10 -extract "$DATA/extract.txt")
[ "$(kept "$out")" = "20" ] && ok "-extract keeps the listed 20 variants" \
                            || bad "-extract keeps the listed 20 variants" "kept $(kept "$out")"

out=$(run preprocess -vcf "$DATA/matched.vcf" -out "$WORK/c_vcf" -ntraits 10 -samples "$DATA/samples.txt")
if grep -q "60 samples" <<<"$out"; then ok "-samples accepts a file of IDs"
else bad "-samples accepts a file of IDs"; fi

# 20 variants x 60 samples x 8 bytes
if [ "$HAVE_PGEN" = "1" ]; then
  run preprocess -vcf "$DATA/matched.vcf" -out "$WORK/d_vcf" -ntraits 10 \
      -extract "$DATA/extract.txt" -samples "$DATA/samples.txt" >/dev/null
  run preprocess -bfile "$DATA/matched" -out "$WORK/d_bed" -ntraits 10 \
      -extract "$DATA/extract.txt" -samples "$DATA/samples.txt" >/dev/null
  size=$(wc -c < "$WORK/d_vcf/chunks/chunk_0.bin" | tr -d ' ')
  if [ "$size" = "9600" ] && cmp -s "$WORK/d_vcf/chunks/chunk_0.bin" "$WORK/d_bed/chunks/chunk_0.bin"; then
    ok "both subsets applied together, readers still agree (9600 bytes)"
  else
    bad "both subsets applied together, readers still agree" "size=$size"
  fi
else
  skip "combined subsetting across readers (built without SAFELD_PGEN)"
fi

# ------------------------------------------------------------ input policy --
echo
echo "input policy"
# 10 records: 3 kept, 2 dup ID, 2 dup locus, 1 all-missing, 1 monomorphic
# (caught by the MAF filter at the default threshold), 1 multiallelic.
out=$(run preprocess -vcf "$DATA/dirty.vcf" -out "$WORK/e" -ntraits 5)
[ "$(kept "$out")" = "3" ] && ok "dirty VCF: exactly 3 of 10 records survive" \
                           || bad "dirty VCF survivors" "kept $(kept "$out")"
grep -q "excluded 4 duplicated" <<<"$out" && ok "duplicates dropped, every copy" \
                                          || bad "duplicates dropped, every copy"
grep -q "excluded 1 non-biallelic" <<<"$out" && ok "multiallelic record skipped" \
                                             || bad "multiallelic record skipped"
grep -q "excluded 1 over the missingness limit" <<<"$out" && ok "all-missing variant reported as missingness" \
                                                          || bad "all-missing variant reported as missingness"
grep -q "excluded 1 below the MAF threshold" <<<"$out" && ok "monomorphic variant reported as MAF exclusion" \
                                                       || bad "monomorphic variant reported as MAF exclusion"

# Every exclusion must be attributed: totals have to balance.
tot=$(grep -oE 'Scanned [0-9,]+' <<<"$out" | tr -d 'Scaned ,')
sum=$(( $(kept "$out") + $(grep -oE 'excluded [0-9,]+' <<<"$out" | tr -d 'excluded ,' | paste -sd+ - | bc) ))
[ "$tot" = "$sum" ] && ok "kept + excluded accounts for every record ($tot)" \
                    || bad "kept + excluded accounts for every record" "scanned $tot, accounted $sum"

# Zero variance is only reachable below the MAF filter, so force -maf 0.
out=$(run preprocess -vcf "$DATA/dirty.vcf" -out "$WORK/e0" -ntraits 5 -maf 0)
grep -q "no variance" <<<"$out" && ok "monomorphic variant dropped by the variance check (-maf 0)" \
                                || bad "monomorphic variant dropped by the variance check (-maf 0)"

# The regression that started all of this: ID-less records must survive.
if grep -qE '^1\s+100\s' "$WORK/e/chunks/chunk_0.meta" && grep -qE '^1\s+200\s' "$WORK/e/chunks/chunk_0.meta"; then
  ok "ID-less ('.') records survive deduplication"
else
  bad "ID-less ('.') records survive deduplication" "they were all dropped, as in the original bug"
fi

echo
echo "dosage field selection"
out=$(run preprocess -vcf "$DATA/sparse_ds.vcf" -out "$WORK/f" -ntraits 5)
grep -q "Dosage source: GT" <<<"$out" && ok "auto reads GT when DS is sparser than GT" \
                                      || bad "auto reads GT when DS is sparser than GT"
out=$(run preprocess -vcf "$DATA/sparse_ds.vcf" -out "$WORK/g" -ntraits 5 -dosage-field DS)
grep -q "filled from" <<<"$out" && ok "forced DS fills absent dosages from hard calls" \
                                || bad "forced DS fills absent dosages from hard calls"
out=$(run preprocess -vcf "$DATA/matched.vcf" -out "$WORK/h" -ntraits 5 -dosage-field DS)
grep -q "declares no DS" <<<"$out" && ok "forced DS on a GT-only file warns and uses GT" \
                                   || bad "forced DS on a GT-only file warns and uses GT"

# --------------------------------------------------------------- pipeline ---
echo
echo "pipeline"
if run simulate -prep "$WORK/a_vcf" -out "$WORK/sim" >/dev/null && \
   run merge -in "$WORK/sim" -out "$WORK/final.vcf" -no-compress >/dev/null; then
  n=$(grep -vc '^#' "$WORK/final.vcf")
  [ "$n" = "60" ] && ok "simulate + merge emit 60 variants" || bad "simulate + merge" "got $n"
  cols=$(grep -m1 '^#CHROM' "$WORK/final.vcf" | awk '{print NF-9}')
  [ "$cols" = "20" ] && ok "output carries 20 trait columns" || bad "trait columns" "got $cols"
else
  bad "simulate + merge run to completion"
fi

# Bad input must fail loudly rather than silently produce nothing.
out=$(run preprocess -vcf "$DATA/matched.vcf" -pfile "$DATA/matched" -out "$WORK/i" 2>&1)
grep -q "exactly one of" <<<"$out" && ok "two inputs at once is rejected" \
                                   || bad "two inputs at once is rejected"

# ------------------------------------------------------------- real data ----
# Skipped unless tests/get_real_testdata.sh has been run. Synthetic fixtures
# cannot produce stale INFO/AF, mixed record types or real allele frequency
# spectra, and every one of those has hidden a bug at some point.
echo
echo "real data (1000 Genomes chr22)"
if [ -f "$REAL/real_ids.vcf" ] && [ -f "$REAL/real.vcf.gz" ]; then
  declare -a names=(vcf bcf)
  run preprocess -vcf "$REAL/real_ids.vcf" -out "$WORK/r_vcf" -ntraits 10 >/dev/null
  run preprocess -vcf "$REAL/real_ids.bcf" -out "$WORK/r_bcf" -ntraits 10 >/dev/null
  if [ "$HAVE_PGEN" = "1" ]; then
    run preprocess -pfile "$REAL/real_ids" -out "$WORK/r_pgen" -ntraits 10 >/dev/null
    run preprocess -bfile "$REAL/real_ids" -out "$WORK/r_bed"  -ntraits 10 >/dev/null
    names+=(pgen bed)
  fi

  n_vcf=$(grep -c '^22' "$WORK/r_vcf/chunks/chunk_0.meta" || echo 0)
  [ "$n_vcf" -gt 100 ] && ok "real VCF keeps $n_vcf variants" \
                       || bad "real VCF keeps a plausible number of variants" "got $n_vcf"

  same=1
  for f in "${names[@]:1}"; do
    cmp -s "$WORK/r_vcf/chunks/chunk_0.bin" "$WORK/r_${f}/chunks/chunk_0.bin" || { same=0; echo "        $f differs"; }
  done
  if [ "$same" = "1" ]; then
    ok "all ${#names[@]} formats agree byte for byte on real data (${names[*]})"
  else
    bad "all formats agree byte for byte on real data"
  fi

  # Every ID in this region is "." — the pathology that broke the original tool,
  # here in real data rather than a constructed fixture.
  out=$(run preprocess -vcf "$REAL/real.vcf.gz" -out "$WORK/r_dots" -ntraits 5)
  k=$(kept "$out")
  [ "${k:-0}" -gt 100 ] && ok "real VCF with no rsIDs at all keeps $k variants" \
                        || bad "real VCF with no rsIDs keeps variants" "kept ${k:-0}"
  grep -q "non-biallelic" <<<"$out" && ok "real multiallelic and symbolic records skipped" \
                                    || bad "real multiallelic and symbolic records skipped"

  # The claim the tool actually makes: LD computed on the synthetic output
  # reproduces LD computed on the real genotypes. Everything else in this file
  # checks that the inputs were read correctly, which is necessary but says
  # nothing about whether the simulation preserves what it is supposed to.
  if command -v plink2 >/dev/null 2>&1; then
    E="$WORK/e2e"; mkdir -p "$E"
    run preprocess -vcf "$REAL/real_ids.vcf" -out "$E/prep" -ntraits 2000 >/dev/null
    run simulate   -prep "$E/prep" -out "$E/sim" -compress >/dev/null
    run merge      -in "$E/sim" -out "$E/safeld.vcf.gz" >/dev/null

    ( cd "$E" && \
      plink2 --vcf safeld.vcf.gz dosage=DS --make-pgen --out SAFELD >/dev/null 2>&1 && \
      plink2 --pfile SAFELD --r2-unphased --ld-window-kb 1000 --ld-window-r2 0 \
             --out SAFELD >/dev/null 2>&1 && \
      plink2 --pfile "$REAL/real_ids" --r2-unphased --ld-window-kb 1000 --ld-window-r2 0 \
             --out ORIGINAL >/dev/null 2>&1 )

    if [ -s "$E/ORIGINAL.vcor" ] && [ -s "$E/SAFELD.vcor" ]; then
      corr=$(awk -f "$HERE/../scripts/ld_cor.awk" "$E/ORIGINAL.vcor" "$E/SAFELD.vcor" 2>/dev/null)
      r=$(grep -oE 'Pearson r  *: [0-9.]+' <<<"$corr" | grep -oE '[0-9.]+$')
      slope=$(grep -oE 'slope \(y~0\+x\) *: [0-9.]+' <<<"$corr" | grep -oE '[0-9.]+$')
      # 2000 traits gives r about 0.998 here; the threshold leaves room for the
      # random trait matrix without admitting a real regression.
      if awk -v r="${r:-0}" 'BEGIN{exit !(r > 0.995)}'; then
        ok "simulated LD reproduces original LD (r=$r, slope=$slope, 2000 traits)"
      else
        bad "simulated LD reproduces original LD" "r=${r:-none}, slope=${slope:-none}"
      fi
      # A slope well below 1 means LD is being lost systematically, which is a
      # different failure from noise and needs saying separately.
      if awk -v s="${slope:-0}" 'BEGIN{exit !(s > 0.97 && s < 1.03)}'; then
        ok "no systematic attenuation (slope $slope within 3% of 1)"
      else
        bad "no systematic attenuation" "slope=${slope:-none}"
      fi
    else
      bad "plink2 produced LD matrices for the end-to-end check"
    fi
  else
    skip "end-to-end LD correlation (plink2 not in PATH)"
  fi

  if [ "$HAVE_PGEN" = "1" ]; then
    a=$(run preprocess -pfile "$REAL/real_ids" -out "$WORK/r_ex_p" -ntraits 5 \
            -extract "$REAL/extract_ids.txt" -samples "$REAL/subset_samples.txt")
    b=$(run preprocess -vcf "$REAL/real_ids.vcf" -out "$WORK/r_ex_v" -ntraits 5 \
            -extract "$REAL/extract_ids.txt" -samples "$REAL/subset_samples.txt")
    if cmp -s "$WORK/r_ex_p/chunks/chunk_0.bin" "$WORK/r_ex_v/chunks/chunk_0.bin"; then
      ok "real data, both subsets, pgen and VCF still agree"
    else
      bad "real data, both subsets, pgen and VCF still agree"
    fi
  fi
else
  skip "real-data checks (run tests/get_real_testdata.sh first)"
fi

echo
printf "  %d passed, %d failed" "$PASS" "$FAIL"
[ "$SKIP" -gt 0 ] && printf ", %d skipped" "$SKIP"
echo
echo
[ "$FAIL" -eq 0 ] || exit 1
