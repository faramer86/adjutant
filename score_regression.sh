#!/bin/sh
# Golden regression harness for adjutant (the --score-only plink2 reduction).
#
# Exercises the full --score / --score-col-nums argument surface and dumps
# every output file into a results directory.
#
# Usage:
#   ./score_regression.sh <adjutant-binary> <output-dir>
#
# Compare two runs with:  diff -r <dir-a> <dir-b>

set -e

BIN=${1:?usage: score_regression.sh <adjutant-binary> <output-dir>}
OUT=${2:?usage: score_regression.sh <adjutant-binary> <output-dir>}

D=/tmp/scoretest
PFILE="Tests/TEST_PHASED_VCF/mt_mixed_phase vzs"

# The test dataset contains duplicate variant IDs and haploid chromosomes.
DUP=ignore-dup-ids
AUTO="--extract $D/autosomal.txt"

rm -rf "$OUT"
mkdir -p "$OUT"

### --- regenerate fixtures (idempotent) ----------------------------------
# Derived from the bundled TEST_PHASED_VCF dataset so the harness is
# self-contained: 2504 samples x 509 variants, incl. multiallelic + haploid.
if [ ! -f "$D/score_multi.txt" ] || [ ! -f "$D/freqs.afreq" ]; then
  mkdir -p "$D"
  cp Tests/TEST_PHASED_VCF/mt_mixed_phase.pvar "$D/mt.pvar"
  # score file: variant ID, ALT1 allele, three deterministic weight columns
  awk 'BEGIN{OFS="\t"; print "ID","A1","W1","W2","W3"}
       !/^#/{split($5,a,","); print $3, a[1], ((NR%17)-8)/3.0, (NR%7)*0.35-1.05, (NR%3)-1}' \
      "$D/mt.pvar" > "$D/score_multi.txt"
  # --q-score-range inputs
  awk 'BEGIN{OFS="\t"} !/^#/{print $3, ((NR*37)%1000)/1000.0}' "$D/mt.pvar" > "$D/qsr_data.txt"
  printf 'lo\t0\t0.2\nmid\t0.2\t0.6\nhi\t0.6\t1.0\nall\t0\t1\n' > "$D/qsr_ranges.txt"
  # --score-list input
  echo "$D/score_multi.txt" > "$D/score_list.txt"
  # companion filter fixtures
  grep -v '^#' Tests/TEST_PHASED_VCF/mt_mixed_phase.psam | awk 'NR%3==0{print $1}' > "$D/keep.txt"
  grep -v '^#' "$D/mt.pvar" | awk 'NR%2==0{print $3}' > "$D/extract.txt"
  # autosomes only: 'dominant'/'recessive'/'variance-standardize' reject haploid
  # chromosomes, and --chr no longer exists in this build.
  grep -v '^#' "$D/mt.pvar" | awk '$1 ~ /^([1-9]|1[0-9]|2[0-2])$/ {print $3}' > "$D/autosomal.txt"
  # --read-freq needs an .afreq.  This used to be generated here with --freq,
  # but --freq is one of the removed commands, so the line failed and `set -e`
  # killed the harness before it printed anything -- exit 8, no diagnostic.
  #
  # It went unnoticed because /tmp/scoretest/freqs.afreq survived on the
  # development machine from before --freq was removed, and the enclosing
  # `if [ ! -f ... ]` then skipped regeneration.  score_check.sh passed there
  # and failed on every clean checkout.  Ship the fixture instead.
  cp Tests/score_fixtures/freqs.afreq "$D/freqs.afreq"
fi

n=0
run() {
  # run <case-name> <extra plink2 args...>
  name=$1
  shift
  n=$((n + 1))
  # shellcheck disable=SC2086
  "$BIN" --pfile $PFILE --out "$OUT/$name" "$@" > "$OUT/$name.stdout" 2>&1 || \
    echo "NONZERO_EXIT $name" >> "$OUT/_failures.txt"
}

### --- core --score behaviour -------------------------------------------
run basic            --score "$D/score_multi.txt" 1 2 3 header $DUP
run cols_all         --score "$D/score_multi.txt" 1 2 3 header $DUP cols=maybefid,fid,maybesid,sid,pheno1,phenos,nallele,denom,dosagesum,scoreavgs,scoresums
run cols_minimal     --score "$D/score_multi.txt" 1 2 3 header $DUP cols=scoresums
run header_read      --score "$D/score_multi.txt" 1 2 3 header-read $DUP
run no_mean_imp      --score "$D/score_multi.txt" 1 2 3 header $DUP no-mean-imputation
run center           --score "$D/score_multi.txt" 1 2 3 header $DUP center
run var_std          --score "$D/score_multi.txt" 1 2 3 header $DUP variance-standardize $AUTO
run dominant         --score "$D/score_multi.txt" 1 2 3 header $DUP dominant $AUTO
run recessive        --score "$D/score_multi.txt" 1 2 3 header $DUP recessive $AUTO
run se_mode          --score "$D/score_multi.txt" 1 2 3 header $DUP se
run list_variants    --score "$D/score_multi.txt" 1 2 3 header $DUP list-variants
run list_variants_zs --score "$D/score_multi.txt" 1 2 3 header $DUP list-variants-zs
run zs_out           --score "$D/score_multi.txt" 1 2 3 header $DUP zs
run implicit_cols    --score "$D/score_multi.txt" header $DUP

### --- multiple score columns -------------------------------------------
run col_nums_range   --score "$D/score_multi.txt" 1 2 header-read $DUP --score-col-nums 3-5
run col_nums_list    --score "$D/score_multi.txt" 1 2 header-read $DUP --score-col-nums 3,5

### --- PLINK1 .bed/.bim/.fam input ---------------------------------------
# Same dataset restricted to biallelic variants (.bim cannot hold multiallelics).
BED="Tests/TEST_PHASED_VCF/mt"
runbed() {
  name=$1; shift
  n=$((n + 1))
  # shellcheck disable=SC2086
  "$BIN" --bfile $BED --out "$OUT/$name" "$@" > "$OUT/$name.stdout" 2>&1 || \
    echo "NONZERO_EXIT $name" >> "$OUT/_failures.txt"
}
runbed bed_basic      $AUTO --score "$D/score_multi.txt" 1 2 3 header $DUP
runbed bed_cols       $AUTO --score "$D/score_multi.txt" 1 2 3 header $DUP cols=maybefid,fid,maybesid,nallele,denom,dosagesum,scoreavgs,scoresums
runbed bed_colnums    $AUTO --score "$D/score_multi.txt" 1 2 header-read $DUP --score-col-nums 3-5
runbed bed_center     $AUTO --score "$D/score_multi.txt" 1 2 3 header $DUP center
runbed bed_keep       $AUTO --score "$D/score_multi.txt" 1 2 3 header $DUP --keep "$D/keep.txt"

### --- interaction with companion filters --------------------------------
run with_threads     --score "$D/score_multi.txt" 1 2 3 header $DUP --threads 1
run with_keep        --score "$D/score_multi.txt" 1 2 3 header $DUP --keep "$D/keep.txt"
run with_extract     --score "$D/score_multi.txt" 1 2 3 header $DUP --extract "$D/extract.txt"
run with_readfreq    --score "$D/score_multi.txt" 1 2 3 header $DUP --read-freq "$D/freqs.afreq" $AUTO

### --- error paths (must still be rejected, with the same message) -------
run err_both_headers --score "$D/score_multi.txt" 1 2 3 header header-read
run err_two_models   --score "$D/score_multi.txt" 1 2 3 header center dominant
run err_dup_ids      --score "$D/score_multi.txt" 1 2 3 header
run err_dom_haploid  --score "$D/score_multi.txt" 1 2 3 header $DUP dominant
run err_colnums_3pos --score "$D/score_multi.txt" 1 2 3 header --score-col-nums 3

### --- flags removed from this build (must be rejected as unrecognised) --
run gone_score_list  --score-list "$D/score_list.txt" 1 2 3 header
run gone_qsr         --score "$D/score_multi.txt" 1 2 3 header $DUP \
                     --q-score-range "$D/qsr_ranges.txt" "$D/qsr_data.txt"
run gone_glm         --glm
run gone_freq        --freq
run gone_export      --export vcf
run gone_pca         --pca
run gone_king        --make-king
run gone_makepgen    --make-pgen
run gone_hardy       --hardy
run gone_vscore      --variant-score "$D/score_multi.txt"

### --- normalise volatile content ---------------------------------------
# Logs and stdout embed timestamps, the hostname, host paths, thread counts
# and timings.  Anything machine-specific must be normalised here, or the
# golden files bake in one developer's environment and every other host reads
# as a --score regression.
for f in "$OUT"/*.log "$OUT"/*.stdout; do
  [ -e "$f" ] || continue
  sed -E -e '/^Start time:/d' \
         -e '/^End time:/d' \
         -e '/threads? \(change/d' \
         -e '/^Logging to /d' \
         -e 's|^Hostname: .*|Hostname: <HOST>|' \
         -e 's/^Random number seed: [0-9]+/Random number seed: <SEED>/' \
         -e 's/[0-9]+ MiB/<MEM> MiB/g' \
         -e 's|/tmp/[A-Za-z0-9_.-]+|<TMP>|g' \
         -e "s|$OUT|<OUT>|g" \
         -e 's/^(PLINK|adjutant) v.*/PLINK <VERSION>/' \
         -e 's/(plink2|adjutant) --help/<PROG> --help/g' \
         -e 's|^Working directory: .*|Working directory: <CWD>|' \
         "$f" > "$f.norm"
  mv "$f.norm" "$f"
done

echo "ran $n cases -> $OUT"
if [ -f "$OUT/_failures.txt" ]; then
  echo "cases with nonzero exit (expected for err_* cases):"
  cat "$OUT/_failures.txt"
fi
