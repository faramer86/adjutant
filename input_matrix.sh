#!/bin/sh
# Input-combination matrix for adjutant.
#
#   ./input_matrix.sh [binary] [test_input-repo]
#
# score_check.sh proves --score is numerically identical to stock plink2.
# This proves the surface around it behaves: that every valid way of naming a
# fileset works, that every invalid one is refused, and that filters compose.
#
# Exit 0 = pass.

set -e

BIN=${1:-./bin/adjutant}
# Default to the vendored dataset, so a fresh clone can run this.  It used to
# point at a directory outside the repo, which meant this suite could not run
# at all on a checkout -- while README told the reader to run it.
T=${2:-Tests/input_matrix}
OUT=${OUT_DIR:-/tmp/adjutant_matrix}

# Fail with a sentence rather than a cascade of awk errors: with the dataset
# absent this used to run on into "can't open file" for every case.
if [ ! -f "$T/pfile/1kg_hg38.pgen" ]; then
  echo "FAIL  no dataset at '$T'" >&2
  echo "      expected $T/{pfile,bfile}/1kg_hg38.* and $T/scores/*.txt" >&2
  echo "      usage: $0 [adjutant-binary] [dataset-dir]" >&2
  exit 1
fi

rm -rf "$OUT"; mkdir -p "$OUT"

P="$T/pfile/1kg_hg38"
B="$T/bfile/1kg_hg38"
W="$T/scores/weights_single.txt"
WM="$T/scores/weights_multi.txt"

pass=0; fail=0

# ok <label> <args...>   -- must succeed and write a .sscore
ok() {
  label=$1; shift
  if "$BIN" "$@" --out "$OUT/$label" > "$OUT/$label.log" 2>&1 && [ -s "$OUT/$label.sscore" ]; then
    pass=$((pass+1)); printf "  ok    %s\n" "$label"
  else
    fail=$((fail+1)); printf "  FAIL  %s\n" "$label"
    grep -E "^Error" "$OUT/$label.log" | head -1 | sed 's/^/          /'
  fi
}

# no <label> <expected-substring> <args...>  -- must be refused, with that message
no() {
  label=$1; want=$2; shift 2
  "$BIN" "$@" --out "$OUT/$label" > "$OUT/$label.log" 2>&1 || true
  if grep -qE "^Error" "$OUT/$label.log" && grep -qF "$want" "$OUT/$label.log"; then
    pass=$((pass+1)); printf "  ok    %s (refused)\n" "$label"
  else
    fail=$((fail+1)); printf "  FAIL  %s: expected refusal containing '%s'\n" "$label" "$want"
    tail -2 "$OUT/$label.log" | sed 's/^/          /'
  fi
}

# okz: same as ok, but the run is expected to write <label>.sscore.zst
okz() {
  label=$1; shift
  if "$BIN" "$@" --out "$OUT/$label" > "$OUT/$label.log" 2>&1 && [ -s "$OUT/$label.sscore.zst" ]; then
    pass=$((pass+1)); printf "  ok    %s\n" "$label"
  else
    fail=$((fail+1)); printf "  FAIL  %s\n" "$label"
  fi
}

# warn_only: run must not error, and must emit the given note/warning
warn_only() {
  label=$1; want=$2; shift 2
  "$BIN" "$@" --out "$OUT/$label" > "$OUT/$label.log" 2>&1 || true
  if grep -qF "$want" "$OUT/$label.log"; then
    pass=$((pass+1)); printf "  ok    %s (declined cleanly)\n" "$label"
  else
    fail=$((fail+1)); printf "  FAIL  %s: expected '%s'\n" "$label" "$want"
  fi
}

echo "--- naming a fileset: complete forms ---"
ok  pfile          --pfile $P            --score $W 1 2 3 header
ok  bfile          --bfile $B            --score $W 1 2 3 header

echo "--- naming a fileset: incomplete (must be refused) ---"
FULL="A full fileset"
no  no_input   "$FULL" --score $W 1 2 3 header
# --pgen and --psam never existed in adjutant
no  gone_pgen  "Unrecognized flag" --pgen $P.pgen --score $W 1 2 3 header
no  gone_psam  "Unrecognized flag" --pfile $P --psam $P.psam --score $W 1 2 3 header
# --bim is now the harmonisation input (--pgs-model); mixing with --score is refused
no  bim_needs_pgs_model "requires --pgs-model" --bim $B.bim --score $W 1 2 3 header

echo "--- conflicting input flags (must be refused) ---"
no  pfile_bfile "conflicts" --pfile $P --bfile $B --score $W 1 2 3 header

echo "--- missing / malformed command ---"
warn_only no_score "No output requested" --pfile $P
no  bad_scorefile "nonexistent.txt does"    --pfile $P --score "$T/nonexistent.txt" 1 2 3 header

echo "--- filters, singly ---"
awk '!/^#/ && NR<=400 {print $1"\t"$2}' $P.psam > "$OUT/keep.txt"
awk '!/^#/ && NR<=300 {print $1}'          $P.psam > "$OUT/keepfam.txt"
awk '!/^#/ && NR<=600 {print $3}'          $P.pvar > "$OUT/extract.txt"
awk '!/^#/ && NR>1800 {print $3}'          $P.pvar > "$OUT/extract2.txt"
ok  f_keep      --pfile $P --keep    "$OUT/keep.txt"    --score $W 1 2 3 header
ok  f_remove    --pfile $P --remove  "$OUT/keep.txt"    --score $W 1 2 3 header
ok  f_extract   --pfile $P --extract "$OUT/extract.txt" --score $W 1 2 3 header
ok  f_exclude   --pfile $P --exclude "$OUT/extract.txt" --score $W 1 2 3 header

echo "--- filters, composed ---"
ok  f_sample_variant --pfile $P --keep "$OUT/keep.txt" --extract "$OUT/extract.txt" --score $W 1 2 3 header
ok  f_extract_exclude --pfile $P --extract "$OUT/extract.txt" --exclude "$OUT/extract2.txt" --score $W 1 2 3 header
ok  f_all_four       --pfile $P --keep "$OUT/keep.txt" --extract "$OUT/extract.txt" \
                     --exclude "$OUT/extract2.txt" --score $W 1 2 3 header
ok  f_on_bfile       --bfile $B --keep "$OUT/keep.txt" --extract "$OUT/extract.txt" --score $W 1 2 3 header

echo "--- --score argument forms ---"
ok  s_explicit    --pfile $P --score $W 1 2 3 header
ok  s_implicit    --pfile $P --score $W header
ok  s_headerread  --pfile $P --score $WM 1 2 header-read --score-col-nums 3-5
ok  s_colnums_csv --pfile $P --score $WM 1 2 header-read --score-col-nums 3,5
ok  s_center      --pfile $P --score $W 1 2 3 header center
ok  s_varstd      --pfile $P --score $W 1 2 3 header variance-standardize
ok  s_dominant    --pfile $P --score $W 1 2 3 header dominant
ok  s_recessive   --pfile $P --score $W 1 2 3 header recessive
ok  s_nomeanimp   --pfile $P --score $W 1 2 3 header no-mean-imputation
ok  s_se          --pfile $P --score $W 1 2 3 header se
ok  s_cols        --pfile $P --score $W 1 2 3 header cols=maybefid,fid,nallele,denom,dosagesum,scoreavgs,scoresums
ok  s_listvars    --pfile $P --score $W 1 2 3 header list-variants
okz s_zs          --pfile $P --score $W 1 2 3 header zs

echo "--- --score modifier conflicts (must be refused) ---"
no  s_two_headers "cannot be used" --pfile $P --score $W 1 2 3 header header-read
no  s_two_models  "mutually exclusive"      --pfile $P --score $W 1 2 3 header center dominant
no  s_colnums_3   "cannot be used when"     --pfile $P --score $W 1 2 3 header --score-col-nums 3

echo "--- scoring companions ---"
ok  c_threads1  --pfile $P --score $W 1 2 3 header --threads 1
ok  c_seed      --pfile $P --score $W 1 2 3 header --seed 1
ok  c_memory    --pfile $P --score $W 1 2 3 header --memory 2048

echo
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
