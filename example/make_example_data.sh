#!/bin/sh
# Rebuild the derived example files from the committed genotypes.
#
#   cd example && ./make_example_data.sh
#
# The genotypes themselves (demo.pgen/.pvar/.psam -- 300 samples, 40 variants)
# were made once with plink2 and are committed, so this script needs only
# adjutant and awk:
#
#   plink2 --dummy 300 40 acgt --make-pgen --seed 20260803 --out demo
#   plink2 --pfile demo --make-bed --out demo    # demo.bed/.bim/.fam
#   plink2 --pfile demo --freq --out demo        # demo.afreq, for --read-freq
#
# Both fileset forms are committed; README uses --bfile.  They score
# identically -- the only difference is that a .fam always carries FID, so
# --bfile output has an FID column and --pfile output does not.
#
# (adjutant has no --freq of its own; --read-freq consumes a file like this.)
#
# Everything else is derived here.  No rand(): the "noise" is a deterministic
# function of the sample index, so this reproduces byte-for-byte on any
# machine, which is what lets README quote real numbers.
set -e
cd "$(dirname "$0")"
A=../bin/adjutant

# ---------------------------------------------------------------------------
# 1. Three candidate PGS models over the 40 variants.
#    PGS_B is the one the phenotype will actually be built from.
# ---------------------------------------------------------------------------
awk 'BEGIN { OFS="\t"; print "ID", "A1", "PGS_A", "PGS_B", "PGS_C" }
     !/^#/ { n++
             print $3, $5, ((n*7)  % 13 - 6) / 10,
                          ((n*3)  % 11 - 5) / 10,
                          ((n*5)  % 7  - 3) / 10 }' demo.pvar > weights.txt

# ---------------------------------------------------------------------------
# 2. Score them, so the phenotype can be built from a real PGS_B.
# ---------------------------------------------------------------------------
$A --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
   --out .tmp_pgs > /dev/null 2>&1

# ---------------------------------------------------------------------------
# 3. HEIGHT is driven by PGS_B; DISEASE is a case/control split on the same
#    signal.  The multipliers (80, 30) are tuned so the effect is realistic:
#    R2_INCR about 0.17 and OR about 2.4 / AUC 0.70.  Larger values make PGS_B
#    separate the classes almost perfectly, and the demo then shows an odds
#    ratio in the millions, which teaches the reader nothing.
# ---------------------------------------------------------------------------
# Two passes: the case/control threshold has to sit at the MEAN of PGS_B, or
# the split is 295/5 and the logistic fit reports a near-separation odds ratio
# in the tens of thousands -- which is a pathology, not a demonstration.
awk -F'\t' '
  NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h)
              if (h == "IID") iid = i; if (h == "PGS_B_AVG") b = i }
            next }
  { sum += $b; n++ }
  END { printf "%.10f\n", sum / n }' .tmp_pgs.sscore > .tmp_mean

awk -F'\t' -v mean="$(cat .tmp_mean)" '
  NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h)
              if (h == "IID") iid = i; if (h == "PGS_B_AVG") b = i }
            print "#IID\tHEIGHT\tDISEASE"; next }
  { n++
    # Deterministic pseudo-noise in [-1, 1), no rand().
    wobble = ((n * 2654435761) % 2000) / 1000 - 1
    height = 170 + 80 * $b + 4 * wobble
    disease = (30 * ($b - mean) + wobble > 0) ? 2 : 1
    printf "%s\t%.2f\t%d\n", $iid, height, disease }' .tmp_pgs.sscore > pheno.txt
rm -f .tmp_mean

# ---------------------------------------------------------------------------
# 4. Covariates, unrelated to the signal.
# ---------------------------------------------------------------------------
awk -F'\t' 'NR == 1 { print "#IID\tAGE\tPC1"; next }
  { n++; printf "%s\t%d\t%.3f\n", $1, 40 + (n * 17) % 40,
           ((n * 31) % 200) / 100 - 1 }' demo.psam > covar.txt

# ---------------------------------------------------------------------------
# 5. Two overlapping variant subsets, for --score-partitions.
# ---------------------------------------------------------------------------
awk 'BEGIN { OFS="\t"; print "SNP", "PATHWAY_A", "PATHWAY_B" }
     !/^#/ { n++; print $3, (n % 2), (n <= 20 ? 1 : 0) }' demo.pvar > partitions.txt

# ---------------------------------------------------------------------------
# 6. Which (phenotype, PGS) pairs are worth fitting, for --pgs-pheno-pairs.
# ---------------------------------------------------------------------------
printf '# phenotype  PGS\nHEIGHT\tPGS_B\nDISEASE\tPGS_B\n' > pairs.txt

# ---------------------------------------------------------------------------
# 7. Harmonisation inputs: a target .bim, and a "published" weight file whose
#    alleles are deliberately not all in the target's orientation -- some
#    swapped, some strand-flipped -- which is the problem --pgs-model solves.
# ---------------------------------------------------------------------------
awk 'BEGIN { OFS="\t" }
     !/^#/ { n++; print $1, $3, 0, $2, $5, $4 }' demo.pvar > target.bim

# Column 6 is the source effect-allele frequency, which --resolve-ambiguous
# needs to decide strand for palindromic (A/T, C/G) variants.  Taken from the
# target's own frequencies and nudged, so the two agree on strand -- the flag's
# whole premise is that source and target come from comparable ancestries.
awk -F'\t' 'FNR == NR { if ($1 !~ /^#/) { split($0, f, "\t"); af[f[2]] = f[5] } next }
     BEGIN { OFS="\t"; comp["A"]="T"; comp["T"]="A"; comp["C"]="G"; comp["G"]="C" }
     FNR == 1 { print "chr", "pos", "effect_allele", "other_allele", "beta", "effect_allele_freq"; next }
     !/^#/ { n++
             ea = $5; oa = $4; frq = af[$3] + 0
             if (n % 3 == 1) { t = ea; ea = oa; oa = t; frq = 1 - frq }   # swapped
             else if (n % 3 == 2) { ea = comp[ea]; oa = comp[oa] }        # flipped
             printf "%s\t%s\t%s\t%s\t%s\t%.4f\n", $1, $2, ea, oa,
                    ((n * 3) % 11 - 5) / 10, frq }' demo.afreq demo.pvar > published_pgs.txt

rm -f .tmp_pgs.sscore .tmp_pgs.log
echo "example data rebuilt:"
ls -1 weights.txt pheno.txt covar.txt partitions.txt pairs.txt target.bim published_pgs.txt | sed 's/^/  /'
