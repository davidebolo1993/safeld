#!/usr/bin/env bash
#
# diagnose_vcf.sh - check whether an input VCF has the properties that make
# safeld and the original R implementation disagree.
#
# The R reference (safe.ld.vcf) was written for Minimac-style imputed dosage
# VCFs: biallelic, DS Number=1 always present, unique CHROM:POS IDs, INFO/AF on
# every record. On that input safeld and the R code select identical variants and
# produce identical dosages. Everything below is a way the input can depart from
# that, and each departure is handled differently by the two implementations.
#
# Usage: ./diagnose_vcf.sh input.vcf.gz [n_records]
#
set -uo pipefail

VCF=${1:?usage: diagnose_vcf.sh input.vcf.gz [n_records]}
N=${2:-200000}

command -v bcftools >/dev/null || { echo "bcftools not found in PATH" >&2; exit 1; }

echo "=============================================================="
echo "  $VCF"
echo "  (first $N records)"
echo "=============================================================="

echo
echo "--- header declarations ------------------------------------"
bcftools view -h "$VCF" 2>/dev/null | grep -E '^##(FORMAT=<ID=(DS|GT)|INFO=<ID=AF)' || \
  echo "  !! no DS/GT FORMAT or AF INFO declaration found"

DS_NUMBER=$(bcftools view -h "$VCF" 2>/dev/null | sed -n 's/.*FORMAT=<ID=DS,Number=\([^,]*\).*/\1/p' | head -1)
if [ -n "$DS_NUMBER" ] && [ "$DS_NUMBER" != "1" ]; then
  echo
  echo "  note: DS is Number=$DS_NUMBER. For Number=A this is one value per ALT"
  echo "        allele, so biallelic records still carry exactly one dosage per"
  echo "        sample and safeld reads them normally. It only matters for"
  echo "        multiallelic records (counted below)."
fi

TMP=$(mktemp); trap 'rm -f "$TMP"' EXIT
bcftools view -H "$VCF" 2>/dev/null | head -"$N" | cut -f1-8 > "$TMP"
TOTAL=$(wc -l < "$TMP")
[ "$TOTAL" -eq 0 ] && { echo "no records read"; exit 1; }

pct() { awk -v a="$1" -v b="$TOTAL" 'BEGIN{printf "%.2f%%", (b?100*a/b:0)}'; }
report() { # label count explanation
  printf "  %-34s %8d  (%s)\n" "$1" "$2" "$(pct "$2")"
  [ "$2" -gt 0 ] && printf "      %s\n" "$3"
}

echo
echo "--- variant properties -------------------------------------"
printf "  %-34s %8d\n" "records examined" "$TOTAL"

MULTI=$(cut -f5 "$TMP" | grep -c ',' || true)
report "multiallelic (ALT has ',')" "$MULTI" \
  "safeld skips these; the R code drops them via a non-numeric INFO/AF. Split with: bcftools norm -m -any"

NOID=$(cut -f3 "$TMP" | grep -cx '\.' || true)
report "ID is '.'" "$NOID" \
  "the ORIGINAL safeld dropped ALL of these (they shared one dedup key). Fixed to key on locus."

DUPID=$(cut -f3 "$TMP" | grep -vx '\.' | sort | uniq -d | wc -l | tr -d ' ')
report "distinct duplicated IDs" "$DUPID" \
  "every copy is dropped, by both implementations"

DUPLOCUS=$(cut -f1,2,4,5 "$TMP" | sort | uniq -d | wc -l | tr -d ' ')
report "distinct duplicated loci" "$DUPLOCUS" \
  "duplicated rows in the genotype matrix if not removed"

NOAF=$(awk -F'\t' '$8 !~ /(^|;)AF=/ {n++} END{print n+0}' "$TMP")
report "no INFO/AF" "$NOAF" \
  "the R code drops these (AF becomes NA); safeld computes AF from the dosages instead"

MULTIAF=$(awk -F'\t' 'match($8,/(^|;)AF=[^;]*/){s=substr($8,RSTART,RLENGTH); if (s ~ /,/) n++} END{print n+0}' "$TMP")
report "INFO/AF has multiple values" "$MULTIAF" \
  "the R code drops these; safeld used only the first ALT's AF"

echo
echo "--- missing genotypes (sampled) ----------------------------"
SAMP=$(( TOTAL < 2000 ? TOTAL : 2000 ))
MISS=$(bcftools view -H "$VCF" 2>/dev/null | head -"$SAMP" | \
  awk -F'\t' '{for(i=10;i<=NF;i++){split($i,a,":"); if(a[1]=="."||a[1]~/^\.[\/|]/||a[1]~/[\/|]\.$/) m++; t++}} END{printf "%d %d", m+0, t+0}')
MCOUNT=${MISS% *}; MTOTAL=${MISS#* }
if [ "${MTOTAL:-0}" -gt 0 ]; then
  awk -v m="$MCOUNT" -v t="$MTOTAL" -v s="$SAMP" 'BEGIN{
    printf "  missing/half-called calls        %8d / %d  (%.3f%%) over %d records\n", m, t, 100*m/t, s
    if (m > 0) {
      print "      Mean-imputation attenuates pairwise r2 roughly in proportion to"
      print "      the call rate, which shows up as distinct lower bands in an"
      print "      r2(safeld) vs r2(original) plot. The R reference has no"
      print "      missingness handling at all, so it cannot reproduce this."
    }
  }'
fi

echo
echo "--- verdict ------------------------------------------------"
if [ "$MULTI" -eq 0 ] && [ "$NOID" -eq 0 ] && [ "$NOAF" -eq 0 ] && [ "${MCOUNT:-0}" -eq 0 ]; then
  echo "  Clean imputed-dosage input: safeld and the R reference should agree."
  echo "  If LD still disagrees, look at the number of traits, not the parsing."
else
  echo "  This input departs from what the R reference assumes (see flags above)."
  echo "  Those records are handled differently by the two implementations, so"
  echo "  compare variant SETS first:"
  echo "    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' safeld_out.vcf.gz | sort > a"
  echo "    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' r_out.vcf       | sort > b"
  echo "    comm -3 a b | head"
fi
echo
