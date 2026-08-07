#!/usr/bin/env bash
#
# get_real_testdata.sh - build a small real-data test set in every input format.
#
# Downloads a 100 kb region of 1000 Genomes phase 3 chromosome 22 and converts
# it to VCF, BCF, .pgen/.pvar/.psam and .bed/.bim/.fam. All four then have to
# produce identical matrices, which is the check that caught a live bug the
# synthetic fixtures could not: plink2 recomputes AC and AN when subsetting
# samples but leaves INFO/AF at its original value, so anything filtering on
# INFO/AF silently uses a frequency the data does not have.
#
# The data is not committed. It is ~600 kb and public (1000 Genomes imposes no
# restrictions on use), but keeping it out of the repository keeps the clone
# small and the provenance explicit.
#
# Requires: bcftools, plink2 and network access.
# Usage: tests/get_real_testdata.sh [OUTDIR]
#
set -euo pipefail

OUT="${1:-$(cd "$(dirname "$0")" && pwd)/realdata}"
REGION="22:17060000-17160000"
URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
N_SAMPLES=200

for tool in bcftools plink2; do
  command -v "$tool" >/dev/null || { echo "$tool not found in PATH" >&2; exit 2; }
done

mkdir -p "$OUT"
cd "$OUT"

if [ ! -s raw.vcf.gz ]; then
  echo "fetching $REGION from 1000 Genomes phase 3 ..."
  bcftools view -r "$REGION" "$URL" -Oz -o raw.vcf.gz
fi

# Subset samples, then drop variants that are monomorphic in the subset. This
# is deliberately the same shape as a real analysis: a cohort subset of a larger
# release, which is exactly the situation where INFO/AF goes stale.
# Not "bcftools query -l ... | head": head closes the pipe, bcftools takes
# SIGPIPE, and pipefail then aborts the script with no output at all.
bcftools query -l raw.vcf.gz > all_samples.txt
head -"$N_SAMPLES" all_samples.txt > keep_samples.txt
bcftools view -S keep_samples.txt raw.vcf.gz -Ou | bcftools view -e 'AC=0' -Oz -o real.vcf.gz
bcftools index -f -t real.vcf.gz

# plink cannot represent multiallelic records in .bed and needs unique IDs, so
# the cross-format set is biallelic with CHROM:POS:REF:ALT identifiers, matching
# what a real pipeline produces. real.vcf.gz keeps its multiallelic records and
# its "." IDs for the pathology tests.
plink2 --vcf real.vcf.gz --max-alleles 2 \
       --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 60 truncate \
       --make-pgen --out real_ids >/dev/null
plink2 --pfile real_ids --make-bed   --out real_ids >/dev/null
plink2 --pfile real_ids --export vcf --out real_ids >/dev/null
plink2 --pfile real_ids --export bcf --out real_ids >/dev/null

# An extract list in the format a real pipeline uses.
grep -v '^#' real_ids.pvar | awk 'NR % 3 == 1 {print $3}' > extract_ids.txt
head -60 keep_samples.txt > subset_samples.txt

cat <<EOF

real test data in $OUT

  real.vcf.gz     $(bcftools view -H real.vcf.gz | wc -l | tr -d ' ') variants, ID column is "." throughout,
                  includes multiallelic and symbolic (<CN2>) records
  real_ids.*      $(grep -vc '^#' real_ids.pvar) biallelic variants x $(grep -vc '^#' real_ids.psam) samples as
                  vcf, bcf, pgen/pvar/psam and bed/bim/fam
  extract_ids.txt $(wc -l < extract_ids.txt | tr -d ' ') IDs
  subset_samples.txt $(wc -l < subset_samples.txt | tr -d ' ') sample IDs

run the suite against it:
  tests/run_tests.sh <safeld> <workdir> $OUT
EOF
