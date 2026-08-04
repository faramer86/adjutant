#!/bin/sh
# Differential test: adjutant --score vs stock plink2, on data neither has seen.
#
# The contract adjutant must hold is "byte-identical to stock plink2", so the
# oracle here is plink2 itself rather than a stored expectation.  Every case
# runs both binaries with identical arguments and compares the .sscore bytes.
#
# Data comes from `plink2 --dummy` with fixed seeds: fresh, but regenerable.
set -e
# Resolve the repo root BEFORE changing directory: $0 may be a relative path,
# and once we cd into the scratch dir it no longer resolves.
ROOT=$(cd "$(dirname "$0")/.." && pwd)
A=${ADJUTANT:-$ROOT/bin/adjutant}

# Generated filesets go to a scratch dir, never into the repo.
WORK=${WORK_DIR:-/tmp/adjutant_score_diff}
mkdir -p "$WORK"
cd "$WORK"
D=$(pwd)
pass=0; fail=0; skip=0
FAILED=""

# ---------------------------------------------------------------------------
# Weight files are built from the dataset's own variant list, with a
# deterministic LCG rather than rand(), so a rerun on another awk gives the
# same coefficients and a diff means a real difference.
# ---------------------------------------------------------------------------
mkweights() {  # mkweights <id-col> <allele-col> <infile> <outfile> <seed>
  awk -v ic="$1" -v ac="$2" -v s="$5" 'BEGIN{OFS="\t"} !/^#/ {
        s=(1103515245*s+12345)%2147483648; w1=(s/2147483648-0.5)*2
        s=(1103515245*s+12345)%2147483648; w2=(s/2147483648-0.5)*2
        s=(1103515245*s+12345)%2147483648; w3=(s/2147483648-0.5)*2
        print $ic,$ac,sprintf("%.6f",w1),sprintf("%.6f",w2),sprintf("%.6f",w3)
      }' "$3" > "$4.body"
  printf 'ID\tA1\tW1\tW2\tW3\n' > "$4"
  cat "$4.body" >> "$4"
  cp "$4.body" "$4.noheader"
  rm -f "$4.body"
}

# ---------------------------------------------------------------------------
# run <label> <args...>
#   Both binaries, identical arguments, byte-compare the scoring output.
#   A case both refuse the same way is reported as skip, not pass: it proves
#   nothing about scoring, and silently counting it as pass is how a matrix
#   ends up green while testing nothing.
# ---------------------------------------------------------------------------
run() {
  label=$1; shift
  ae=0; pe=0
  "$A"   "$@" --out "$D/adj_$label" > "$D/adj_$label.out" 2>&1 || ae=$?
  plink2 "$@" --out "$D/plk_$label" > "$D/plk_$label.out" 2>&1 || pe=$?

  if [ "$ae" -ne 0 ] || [ "$pe" -ne 0 ]; then
    if [ "$ae" -ne 0 ] && [ "$pe" -ne 0 ]; then
      skip=$((skip+1)); printf "  skip  %-34s both refused\n" "$label"
    else
      fail=$((fail+1)); FAILED="$FAILED $label"
      printf "  FAIL  %-34s exit differs (adjutant=%s plink2=%s)\n" "$label" "$ae" "$pe"
      grep -m1 -E "^Error" "$D/adj_$label.out" "$D/plk_$label.out" 2>/dev/null | sed 's/^/          /'
    fi
    return
  fi

  # Locate the scoring output: plain, or zstd-compressed under 'zs'.
  if [ -f "$D/adj_$label.sscore" ]; then
    af="$D/adj_$label.sscore"; pf="$D/plk_$label.sscore"
  elif [ -f "$D/adj_$label.sscore.zst" ]; then
    # Compare decompressed bytes.  The two binaries bundle different zstd
    # builds, so identical content can compress to different files; that is
    # not a scoring difference and must not be reported as one.
    zstd -dcq "$D/adj_$label.sscore.zst" > "$D/adj_$label.unz" 2>/dev/null || { skip=$((skip+1)); printf "  skip  %-34s no zstd CLI\n" "$label"; return; }
    zstd -dcq "$D/plk_$label.sscore.zst" > "$D/plk_$label.unz" 2>/dev/null
    af="$D/adj_$label.unz"; pf="$D/plk_$label.unz"
  else
    fail=$((fail+1)); FAILED="$FAILED $label"
    printf "  FAIL  %-34s no .sscore produced\n" "$label"; return
  fi

  if ! cmp -s "$af" "$pf"; then
    fail=$((fail+1)); FAILED="$FAILED $label"
    printf "  FAIL  %-34s .sscore differs\n" "$label"
    diff "$af" "$pf" | head -4 | sed 's/^/          /'
    return
  fi

  # list-variants writes a second file; it is part of the output contract.
  for ext in sscore.vars; do
    if [ -f "$D/adj_$label.$ext" ] && ! cmp -s "$D/adj_$label.$ext" "$D/plk_$label.$ext"; then
      fail=$((fail+1)); FAILED="$FAILED $label($ext)"
      printf "  FAIL  %-34s .%s differs\n" "$label" "$ext"; return
    fi
  done
  pass=$((pass+1)); printf "  ok    %-34s\n" "$label"
}

# ---------------------------------------------------------------------------
# Fresh datasets.  Fixed seeds, so a failure is reproducible.
# ---------------------------------------------------------------------------
echo "Generating fresh datasets with plink2 --dummy ..."
gen() { plink2 --seed "$1" --dummy $2 --out "$3" $4 > /dev/null 2>&1; }
gen 20260803 "400 800 acgt scalar-pheno pheno-ct=2"        dsA     --make-pgen
gen 20260803 "400 800 acgt scalar-pheno pheno-ct=2"        dsAbed  --make-bed
gen 20260804 "300 500 0.05 acgt scalar-pheno"              dsB     --make-pgen
gen 20260805 "250 400 acgt scalar-pheno dosage-freq=0.4"   dsC     --make-pgen
gen 20260806 "40 50 acgt scalar-pheno"                     dsD     --make-pgen
gen 20260807 "200 300 acgt"                                dsE     --make-pgen
gen 20260808 "200 300 acgt scalar-pheno phase-freq=0.3"    dsF     --make-pgen

for d in dsA dsB dsC dsD dsE dsF; do
  mkweights 3 5 "$d.pvar" "w_$d.txt" 999 2>/dev/null || true
done
mkweights 2 5 dsAbed.bim w_dsAbed.txt 999

# --read-freq needs a real .afreq, produced by the oracle, not by us.
plink2 --pfile dsA --freq --out dsA > /dev/null 2>&1

# Filter lists: first/last halves, so each filter actually removes something.
awk 'NR>1 && NR<=120 {print $1}'  dsA.psam > keep.txt
awk 'NR>1 && NR<=60  {print $1}'  dsA.psam > remove.txt
awk '!/^#/ && NR<=300 {print $3}' dsA.pvar > extract.txt
awk '!/^#/ && NR<=150 {print $3}' dsA.pvar > exclude.txt

echo "  $(ls ds*.pgen ds*.bed 2>/dev/null | wc -l | tr -d ' ') filesets, $(ls w_ds*.txt | wc -l | tr -d ' ') weight files"
echo

echo "--- input routes, both fileset forms ---"
run pfile_basic  --pfile dsA    --score w_dsA.txt    1 2 header-read --score-col-nums 3
run bfile_basic  --bfile dsAbed --score w_dsAbed.txt 1 2 header-read --score-col-nums 3

echo "--- header handling ---"
run hdr_read     --pfile dsA --score w_dsA.txt          1 2 header-read --score-col-nums 3
run hdr_skip     --pfile dsA --score w_dsA.txt          1 2 header      --score-col-nums 3
run hdr_none     --pfile dsA --score w_dsA.txt.noheader 1 2             --score-col-nums 3

echo "--- coefficient column selection ---"
run cols_single  --pfile dsA --score w_dsA.txt 1 2 header-read --score-col-nums 3
run cols_two     --pfile dsA --score w_dsA.txt 1 2 header-read --score-col-nums 3-4
run cols_three   --pfile dsA --score w_dsA.txt 1 2 header-read --score-col-nums 3-5
run cols_list    --pfile dsA --score w_dsA.txt 1 2 header-read --score-col-nums 3,5
run cols_default --pfile dsA --score w_dsA.txt 1 2 header-read

echo "--- genotype transforms ---"
run xf_center    --pfile dsA --score w_dsA.txt 1 2 header-read center               --score-col-nums 3
run xf_varstd    --pfile dsA --score w_dsA.txt 1 2 header-read variance-standardize --score-col-nums 3
run xf_dominant  --pfile dsA --score w_dsA.txt 1 2 header-read dominant             --score-col-nums 3
run xf_recessive --pfile dsA --score w_dsA.txt 1 2 header-read recessive            --score-col-nums 3

echo "--- missing-data handling (dsB has 5% missing genotypes) ---"
run miss_default --pfile dsB --score w_dsB.txt 1 2 header-read                    --score-col-nums 3
run miss_nomean  --pfile dsB --score w_dsB.txt 1 2 header-read no-mean-imputation --score-col-nums 3
run miss_center  --pfile dsB --score w_dsB.txt 1 2 header-read center             --score-col-nums 3
run miss_varstd  --pfile dsB --score w_dsB.txt 1 2 header-read variance-standardize --score-col-nums 3
run miss_readfrq --pfile dsA --score w_dsA.txt 1 2 header-read --read-freq dsA.afreq --score-col-nums 3

echo "--- dosages (dsC) and phase (dsF) ---"
run dos_default  --pfile dsC --score w_dsC.txt 1 2 header-read                      --score-col-nums 3
run dos_center   --pfile dsC --score w_dsC.txt 1 2 header-read center               --score-col-nums 3
run dos_varstd   --pfile dsC --score w_dsC.txt 1 2 header-read variance-standardize --score-col-nums 3
run dos_nomean   --pfile dsC --score w_dsC.txt 1 2 header-read no-mean-imputation   --score-col-nums 3
run dos_dominant --pfile dsC --score w_dsC.txt 1 2 header-read dominant             --score-col-nums 3
run phase_basic  --pfile dsF --score w_dsF.txt 1 2 header-read                      --score-col-nums 3

echo "--- standard errors, variant listing, compression ---"
run se_mode      --pfile dsA --score w_dsA.txt 1 2 header-read se             --score-col-nums 3
run listvars     --pfile dsA --score w_dsA.txt 1 2 header-read list-variants  --score-col-nums 3
run zs_out       --pfile dsA --score w_dsA.txt 1 2 header-read zs             --score-col-nums 3

echo "--- column sets ---"
run cs_scoresums --pfile dsA --score w_dsA.txt 1 2 header-read cols=+scoresums --score-col-nums 3
run cs_denom     --pfile dsA --score w_dsA.txt 1 2 header-read cols=+denom     --score-col-nums 3
run cs_minimal   --pfile dsA --score w_dsA.txt 1 2 header-read cols=scoreavgs  --score-col-nums 3
run cs_fidsid    --pfile dsA --score w_dsA.txt 1 2 header-read cols=fid,sid,scoreavgs --score-col-nums 3
run cs_pheno1    --pfile dsA --score w_dsA.txt 1 2 header-read cols=pheno1,scoreavgs  --score-col-nums 3
run cs_phenos    --pfile dsA --score w_dsA.txt 1 2 header-read cols=phenos,scoreavgs  --score-col-nums 3
run cs_nallele   --pfile dsA --score w_dsA.txt 1 2 header-read cols=nallele,dosagesum,scoreavgs --score-col-nums 3

echo "--- filters composed with scoring ---"
run flt_keep     --pfile dsA --keep keep.txt       --score w_dsA.txt 1 2 header-read --score-col-nums 3
run flt_remove   --pfile dsA --remove remove.txt   --score w_dsA.txt 1 2 header-read --score-col-nums 3
run flt_extract  --pfile dsA --extract extract.txt --score w_dsA.txt 1 2 header-read --score-col-nums 3
run flt_exclude  --pfile dsA --exclude exclude.txt --score w_dsA.txt 1 2 header-read --score-col-nums 3
run flt_both     --pfile dsA --keep keep.txt --extract extract.txt --score w_dsA.txt 1 2 header-read --score-col-nums 3
run flt_bed_keep --bfile dsAbed --keep keep.txt    --score w_dsAbed.txt 1 2 header-read --score-col-nums 3

echo "--- edge sizes and phenotype shapes ---"
run tiny_ds      --pfile dsD --score w_dsD.txt 1 2 header-read --score-col-nums 3
run tiny_varstd  --pfile dsD --score w_dsD.txt 1 2 header-read variance-standardize --score-col-nums 3
run binpheno     --pfile dsE --score w_dsE.txt 1 2 header-read --score-col-nums 3
run binpheno_ph  --pfile dsE --score w_dsE.txt 1 2 header-read cols=phenos,scoreavgs --score-col-nums 3

echo "--- combinations that stress several paths at once ---"
run mix_1 --pfile dsB --score w_dsB.txt 1 2 header-read center no-mean-imputation --score-col-nums 3-5
run mix_2 --pfile dsC --score w_dsC.txt 1 2 header-read variance-standardize se    --score-col-nums 3-4
run mix_3 --pfile dsA --keep keep.txt --extract extract.txt --score w_dsA.txt 1 2 header-read center cols=+scoresums --score-col-nums 3-5
run mix_4 --pfile dsA --score w_dsA.txt 1 2 header-read list-variants cols=+denom --read-freq dsA.afreq --score-col-nums 3-4
run mix_5 --bfile dsAbed --remove remove.txt --score w_dsAbed.txt 1 2 header-read dominant --score-col-nums 3-5

echo
echo "======================================================================"
echo "  --score differential vs stock plink2:  $pass ok, $fail FAIL, $skip skipped"
[ -n "$FAILED" ] && echo "  failing:$FAILED"
echo "======================================================================"
[ "$fail" -eq 0 ]
