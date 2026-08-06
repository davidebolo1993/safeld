#!/usr/bin/env bash
#
# ld_gap_stats.sh - aggregate statistics explaining why safeld LD can sit below
# the original LD on a given input, when the input has no missingness, no
# multiallelic records and no duplicate IDs.
#
# WHAT IT PRINTS
#   Only aggregate numbers: counts, means, and histograms over variants. No
#   sample identifiers, no per-sample genotypes, no per-variant records. The
#   output is a few dozen lines and is safe to copy off a secured host by hand.
#
# WHAT IT TESTS
#   safeld consumes the DS (dosage) field. If the "original" LD it is being
#   compared against was computed from GT hard calls, then any variant whose
#   dosage is shrunk relative to its hard call - which is what imputation does to
#   poorly imputed sites - loses LD in safeld while keeping it in the original.
#   The attenuation is predictable: for a pair of variants i,j,
#
#       r2_dosage(i,j)  ~=  r2_hardcall(i,j) * rho2(i) * rho2(j)
#
#   where rho2(v) = cor(DS_v, GT_v)^2. So the DISTRIBUTION of rho2 below predicts
#   where the bands in an r2-vs-r2 plot land. If rho2 is bimodal (array-genotyped
#   sites at ~1.0, imputed sites lower), you get discrete bands. If rho2 is ~1.0
#   everywhere, this explanation is dead and the cause is elsewhere.
#
# REQUIREMENTS: bcftools and awk. Nothing else.
#
# USAGE
#   ./ld_gap_stats.sh INPUT.vcf[.gz] [N_RECORDS] [REGION]
#
#   N_RECORDS  variants to sample (default 5000; use 0 for all)
#   REGION     optional bcftools region, e.g. chr1:1-5000000
#
set -uo pipefail

VCF=${1:?usage: ld_gap_stats.sh INPUT.vcf[.gz] [N_RECORDS] [REGION]}
NREC=${2:-5000}
REGION=${3:-}

command -v bcftools >/dev/null || { echo "bcftools not found in PATH" >&2; exit 1; }

TMPQ=$(mktemp); trap 'rm -f "$TMPQ"' EXIT

echo "=================================================================="
echo " safeld LD-gap statistics"
echo " file    : $(basename "$VCF")"
echo " sampled : ${NREC:-all} records${REGION:+  region: $REGION}"
echo " date    : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "=================================================================="

# ---------------------------------------------------------------- header ----
echo
echo "[1] HEADER"
NSAMP=$(bcftools query -l "$VCF" 2>/dev/null | wc -l | tr -d ' ')
echo "    samples in file          : $NSAMP"
for tag in DS GT GP HDS; do
  line=$(bcftools view -h "$VCF" 2>/dev/null | grep -m1 "FORMAT=<ID=$tag," || true)
  [ -n "$line" ] && echo "    FORMAT $tag               : present ($(echo "$line" | sed -n 's/.*Number=\([^,]*\).*/Number=\1/p'))"
done
INFO_IDS=$(bcftools view -h "$VCF" 2>/dev/null | sed -n 's/^##INFO=<ID=\([^,]*\).*/\1/p' | tr '\n' ' ')
echo "    INFO fields declared     : ${INFO_IDS:-<none>}"
echo "    -> an imputation-quality field (R2/INFO/ER2) is what the chr22 test"
echo "       data was filtered on (R2>=0.8). Its absence here means no such"
echo "       filter was applied."

# ---------------------------------------------------------------- stream ----
# No arrays here: `set -u` plus an empty array is an error on bash 3.x, which
# would silently yield an all-zero report.
STREAM() {
  if [ -n "$REGION" ]; then
    bcftools query -r "$REGION" -f '[\t%GT=%DS]\n' "$VCF" 2>/dev/null
  else
    bcftools query -f '[\t%GT=%DS]\n' "$VCF" 2>/dev/null
  fi
}

if [ "$NREC" = "0" ]; then
  STREAM > "$TMPQ"
else
  STREAM 2>/dev/null | head -"$NREC" > "$TMPQ"
fi

if [ ! -s "$TMPQ" ]; then
  echo
  echo "  !! bcftools produced no output. Check that the file has both GT and DS," >&2
  echo "     and that the region (if given) exists. Aborting rather than" >&2
  echo "     reporting zeros." >&2
  exit 1
fi

awk '
function bin(x, edges, n,   i) { for (i = 1; i <= n; i++) if (x < edges[i]) return i; return n + 1 }
BEGIN {
  FS = "\t"
  ne = 7
  split("0.50 0.70 0.80 0.90 0.95 0.99 0.999", E, " ")
  split("0.50 0.70 0.80 0.90 0.95 0.99 0.999", V, " ")
}
{
  n = 0; sh = 0; sh2 = 0; sd = 0; sd2 = 0; shd = 0
  multi_ds = 0
  for (i = 2; i <= NF; i++) {
    p = index($i, "=")
    if (p == 0) continue
    gt = substr($i, 1, p - 1)
    ds = substr($i, p + 1)

    if (index(ds, ",")) { multi_ds = 1; continue }     # Number=A on a multiallelic
    if (ds == "." || ds == "") continue
    if (index(gt, ".")) continue                        # missing/half hard call

    # hard-call ALT dosage
    nalt = 0; np = 0; a = gt
    gsub(/\|/, "/", a)
    m = split(a, AL, "/")
    for (k = 1; k <= m; k++) { if (AL[k] == ".") { np = -1; break } ; np++; if (AL[k] + 0 > 0) nalt++ }
    if (np <= 0) continue
    h = 2.0 * nalt / np
    d = ds + 0

    n++; sh += h; sh2 += h * h; sd += d; sd2 += d * d; shd += h * d
    total_ds++
    if (d - int(d + 0.5) < 0.01 && int(d + 0.5) - d < 0.01) near_int++
  }
  if (multi_ds) { nmulti++ }
  if (n < 10) { nskip++; next }

  nvar++
  mh = sh / n; md = sd / n
  vh = sh2 / n - mh * mh
  vd = sd2 / n - md * md
  cv = shd / n - mh * md

  af = md / 2.0
  maf = (af < 0.5 ? af : 1 - af)
  sum_maf += maf
  if (maf < 0.01) nmaf01++
  else if (maf < 0.05) nmaf05++

  if (vh <= 0 || vd <= 0) { ndeg++; next }

  rho2 = (cv * cv) / (vh * vd)
  if (rho2 > 1) rho2 = 1
  ratio = vd / vh

  sum_rho2 += rho2; nrho++
  if (rho2 >= 0.999) n_perfect++
  RB[bin(rho2, E, ne)]++
  VB[bin(ratio, V, ne)]++
}
END {
  printf "\n[2] VARIANTS\n"
  printf "    variants usable          : %d\n", nvar
  printf "    skipped (<10 calls)      : %d\n", nskip + 0
  printf "    zero-variance (GT or DS) : %d   <- safeld now drops these\n", ndeg + 0
  printf "    records with comma DS    : %d   <- Number=A on multiallelic\n", nmulti + 0

  printf "\n[3] DOSAGE vs HARD CALL   rho2 = cor(DS, GT)^2 per variant\n"
  if (nrho == 0) { printf "    no usable variants\n"; exit }
  printf "    mean rho2                : %.4f\n", sum_rho2 / nrho
  printf "    variants with rho2>0.999 : %d  (%.1f%%)  <- effectively hard calls\n", n_perfect + 0, 100 * n_perfect / nrho
  printf "    distribution:\n"
  lab[1] = "rho2 <  0.50"; lab[2] = "0.50 - 0.70"; lab[3] = "0.70 - 0.80"
  lab[4] = "0.80 - 0.90"; lab[5] = "0.90 - 0.95"; lab[6] = "0.95 - 0.99"
  lab[7] = "0.99 - 0.999"; lab[8] = "rho2 >= 0.999"
  for (i = 1; i <= 8; i++) {
    c = RB[i] + 0; pc = 100 * c / nrho
    bar = ""; nb = int(pc / 2 + 0.5); for (j = 0; j < nb; j++) bar = bar "#"
    printf "      %-14s %7d  %5.1f%%  %s\n", lab[i], c, pc, bar
  }

  printf "\n[4] VARIANCE RATIO   var(DS)/var(GT) per variant\n"
  vlab[1] = "< 0.50"; vlab[2] = "0.50 - 0.70"; vlab[3] = "0.70 - 0.80"
  vlab[4] = "0.80 - 0.90"; vlab[5] = "0.90 - 0.95"; vlab[6] = "0.95 - 0.99"
  vlab[7] = "0.99 - 0.999"; vlab[8] = ">= 0.999"
  for (i = 1; i <= 8; i++) {
    c = VB[i] + 0; pc = 100 * c / nrho
    printf "      %-14s %7d  %5.1f%%\n", vlab[i], c, pc
  }

  printf "\n[5] DOSAGE GRANULARITY\n"
  printf "    DS values within 0.01 of an integer : %.2f%%\n", (total_ds ? 100 * near_int / total_ds : 0)
  printf "    -> near 100%% means DS is just the hard call and this whole\n"
  printf "       explanation is ruled out. Well below 100%% means genuine\n"
  printf "       imputed dosages, which attenuate LD relative to hard calls.\n"

  printf "\n[6] ALLELE FREQUENCY (from DS)\n"
  printf "    mean MAF                 : %.4f\n", (nvar ? sum_maf / nvar : 0)
  printf "    MAF < 0.01               : %d  (%.1f%%)\n", nmaf01 + 0, 100 * nmaf01 / nvar
  printf "    MAF 0.01 - 0.05          : %d  (%.1f%%)\n", nmaf05 + 0, 100 * nmaf05 / nvar

  printf "\n[7] PREDICTED LD ATTENUATION\n"
  mr = sum_rho2 / nrho
  printf "    If the reference LD came from GT hard calls and safeld from DS,\n"
  printf "    a typical pair loses a factor rho2(i)*rho2(j):\n"
  printf "      typical pair   : %.3f  (r2=0.80 reads as %.3f)\n", mr * mr, 0.80 * mr * mr
  printf "      both perfect   : 1.000  (r2=0.80 reads as 0.800)\n"
  printf "    A gap between these two lines is exactly the banding seen in an\n"
  printf "    r2(safeld) vs r2(original) scatter.\n"

  printf "\n[SUMMARY LINE - copy this back]\n"
  printf "  nvar=%d mean_rho2=%.4f pct_rho2_ge_0.999=%.1f ", nvar, mr, 100 * n_perfect / nrho
  printf "pct_near_int_DS=%.2f mean_maf=%.4f pct_maf_lt_0.01=%.1f degenerate=%d\n", (total_ds ? 100 * near_int / total_ds : 0), (nvar ? sum_maf / nvar : 0), 100 * nmaf01 / nvar, ndeg + 0
}
' "$TMPQ"

echo
echo "=================================================================="
echo " HOW TO READ THIS"
echo "  * [3] mean rho2 near 1.00 and [5] near 100%  -> DS == hard calls."
echo "    The DS/GT explanation is dead; report ntraits used for the run."
echo "  * [3] spread out, or bimodal with a spike at >=0.999 -> genotyped"
echo "    and imputed variants form separate classes. That is the banding,"
echo "    and it is a property of the data, not of the C++ code."
echo "  * [2] zero-variance > 0 -> those variants used to be emitted as a"
echo "    constant dosage of 1.0 for every trait; they are now dropped."
echo "=================================================================="
