#!/bin/sh
# Unit tests for the PGS evaluation numerics, plus an end-to-end CLI check.
#
# Reference values in Tests/pgs_eval_fixtures/ were generated once and
# committed: coefficients and R^2 by statsmodels (generate.py,
# generate_logistic.py), confidence intervals by R (generate_ci.R).  Neither
# Python nor R is needed to run this.
set -e
cd "$(dirname "$0")"

# ---------------------------------------------------------------------------
# Columns are addressed by NAME, never by literal index.
#
# The literal awk indices this replaces were interchangeable with each other:
# in any one .best_pgs most columns are constant down the file or NA on every
# row, and the CI columns substitute for one another -- so "identical on all
# rows", "sorted descending" and "is NA" all pass just as happily on the wrong
# column.  Mutating 16 indices across 12 assertions at once left the suite
# reporting 51 PASS, 0 FAIL.
#
# Reading the index out of the file's own header closes that: a layout change
# moves the name and the index together, and a name that is not in the header
# aborts the run instead of quietly yielding $"".  Case 0 still pins the header
# string itself, which is what makes the names mean anything.
# ---------------------------------------------------------------------------
# saidflat <file> <pattern> -- grep a message with the line wrapping removed.
#
# Every message assertion in this file goes through this, and it is not
# defensive programming -- it is a bug that was shipped and then observed.
#
# adjutant emits through logerrprintfww, which word-wraps to a fixed width, so
# whether a phrase survives on one line depends on how long the interpolated
# path was.  $WDIR comes from `mktemp -d`, whose length differs between
# machines: the "matched no variants" assertion passes under macOS's long
# /var/folders/... prefix and FAILS under a container's short /tmp/tmp.XXXX,
# because the shorter path leaves the line filling up exactly at that phrase.
# A knife-edge, not a platform difference -- which makes any plain grep -q on a
# wrapped message a latent flake that reports a correct binary as broken.
#
# Flattening first tests what the message SAYS rather than where the wrapper
# happened to break it.
saidflat() {
  tr '\n' ' ' < "$1" | tr -s ' ' | grep -q "$2"
}

# need <file> <what was expected> -- turn a missing output into a named failure.
#
# This exists because of a specific way the suite can lie.  Under `set -e` a
# refused adjutant run aborts the script, and so does the `VAR=$(head -1 ...)`
# that reads its output -- both before any assertion prints.  Red-checking the
# --give-corr-matrix guards reported "no failure" for two rounds while the
# mutation was working exactly as intended and the suite was quietly dying at
# line one of the block.  Every run that can fail is followed by a need().
#
# Defined up here with saidflat(), not beside its first heavy user: a helper
# defined halfway down a 3000-line script is a "command not found" waiting for
# whoever adds a case above it.
need() {
  [ -s "$1" ] || { echo "FAIL  $2 (no $1)"; exit 1; }
}

col() {  # col <file> <COLUMN NAME> -> 1-based index
  _idx=$(awk -F'\t' -v want="$2" \
    'NR==1{for(i=1;i<=NF;i++){n=$i; sub(/^#/,"",n); if(n==want){print i; exit}}}' \
    "$1")
  [ -n "$_idx" ] || { echo "FAIL  $1 has no column named '$2'" >&2; exit 1; }
  printf '%s\n' "$_idx"
}

# BLASFLAGS64 is left to the Makefile by default -- its Darwin block sets
# -framework Accelerate, and passing the variable on the command line would
# override that.  Export BLASFLAGS64 to build elsewhere, e.g. on Debian:
#   BLASFLAGS64="-lopenblas -llapacke" ./pgs_eval_check.sh
#
# `adjutant` is built here, not just `pgs_eval_test`.  Most of the cases below
# invoke bin/adjutant, so building only the unit-test binary let the suite
# report green against whatever adjutant happened to be on disk -- including one
# built from different source.  SKIP_BUILD=1 skips only the adjutant build, for
# callers that just built it themselves (tools/verify_and_commit.sh);
# pgs_eval_test is always built because nothing else builds it.
TARGETS="pgs_eval_test"
[ "${SKIP_BUILD:-0}" = "1" ] || TARGETS="$TARGETS adjutant"
make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
     ${BLASFLAGS64:+BLASFLAGS64="$BLASFLAGS64"} \
     LDFLAGS="-ldl -lpthread -lz" $TARGETS > /tmp/pgs_eval_build.log 2>&1 || {
  echo "FAIL  build"; grep -E " error:|ld returned" /tmp/pgs_eval_build.log | head -5; exit 1;
}

if ./bin/pgs_eval_test; then
  echo "PASS  pgs_eval numerics"
else
  echo "FAIL  pgs_eval numerics"; exit 1
fi

# ---------------------------------------------------------------------------
# End-to-end CLI check: --give-best-pgs
# ---------------------------------------------------------------------------
WDIR=$(mktemp -d)
trap 'rm -rf "$WDIR"' EXIT

PSAM="Tests/TEST_PHASED_VCF/mt_mixed_phase.psam"
PFILE_PFX="Tests/TEST_PHASED_VCF/mt_mixed_phase"

awk 'NR==1{print "#IID\tQT1"} NR>1{print $1"\t"(NR%97)/13.0}' \
  "$PSAM" > "$WDIR/pheno.txt"
awk 'NR==1{print "#IID\tC1\tC2"} NR>1{print $1"\t"(NR%7)/3.0"\t"(NR%5)/2.0}' \
  "$PSAM" > "$WDIR/covar.txt"
awk 'BEGIN{print "ID\tA1\tW1\tW2\tW3"} !/^#/{print $3"\t"$5"\t"(NR%7)/10"\t"(NR%5)/10"\t"(NR%3)/10}' \
  Tests/TEST_PHASED_VCF/mt_mixed_phase.pvar | head -400 > "$WDIR/w3.txt"

bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --out "$WDIR/gb" > /dev/null 2>&1

BEST="$WDIR/gb.best_pgs"

# Case 0: the header is exactly this, names and order.
#
# Every other case below resolves columns by name through col(), so this is the
# one place the names and their order are stated rather than looked up.  Without
# it the layout is unchecked -- renaming or reordering any column keeps every
# other case green.  Verified by mutation: renaming #PGS_ID to #MUTANT passed
# 48/48 before this case existed.  When the layout changes, this string is the
# only thing that has to change with it.
EXPECTED_HEADER='#PGS_ID	PHENO	TYPE	N	N_CASE	FOLDS	BETA	OR	SE	P	R2_NULL	R2_FULL	R2_INCR	R2_INCR_L95	R2_INCR_U95	R2L_INCR	R2L_INCR_L95	R2L_INCR_U95	AUC_NULL	AUC_FULL	AUC_INCR'
ACTUAL_HEADER=$(head -1 "$BEST")
if [ "$ACTUAL_HEADER" = "$EXPECTED_HEADER" ]; then
  echo "PASS  give-best-pgs: .best_pgs header exact (names and order)"
else
  echo "FAIL  give-best-pgs: .best_pgs header mismatch"
  echo "  expected: $EXPECTED_HEADER"
  echo "  got:      $ACTUAL_HEADER"; exit 1
fi

# Bind the columns the cases below read.  One binding serves every .best_pgs in
# this file: they all come out of the same writer, and Case 0b asserts the field
# count on the one file that exercises its other (unconverged) branch.
C_PGS_ID=$(col "$BEST" PGS_ID)
C_PHENO=$(col "$BEST" PHENO)
C_TYPE=$(col "$BEST" TYPE)
C_N=$(col "$BEST" N)
C_N_CASE=$(col "$BEST" N_CASE)
C_BETA=$(col "$BEST" BETA)
C_OR=$(col "$BEST" OR)
C_P=$(col "$BEST" P)
C_R2_NULL=$(col "$BEST" R2_NULL)
C_R2_FULL=$(col "$BEST" R2_FULL)
C_R2_INCR=$(col "$BEST" R2_INCR)
C_R2_INCR_L95=$(col "$BEST" R2_INCR_L95)
C_R2_INCR_U95=$(col "$BEST" R2_INCR_U95)
C_R2L_INCR=$(col "$BEST" R2L_INCR)
C_R2L_INCR_L95=$(col "$BEST" R2L_INCR_L95)
C_R2L_INCR_U95=$(col "$BEST" R2L_INCR_U95)
C_AUC_NULL=$(col "$BEST" AUC_NULL)
C_AUC_FULL=$(col "$BEST" AUC_FULL)
C_AUC_INCR=$(col "$BEST" AUC_INCR)

# ---------------------------------------------------------------------------
# Case 0b: every row carries as many fields as the header names, unconverged
# rows included.  The other half of the layout guard: Case 0 pins the names,
# this pins the count, and an awk index only means anything if both hold.
#
# W4 has all-zero weights, so its score is constant, PgsEvalOlsOne rejects the
# column, and the row comes out unconverged -- the one path that writes its
# fields as a fixed run of NA rather than one at a time.  That run was one field
# short for the entire life of the 16-column layout and no case here produced an
# unconverged row to notice.
# ---------------------------------------------------------------------------
awk 'BEGIN{print "ID\tA1\tW1\tW2\tW3\tW4"}
     !/^#/{print $3"\t"$5"\t"(NR%7)/10"\t"(NR%5)/10"\t"(NR%3)/10"\t0"}' \
  Tests/TEST_PHASED_VCF/mt_mixed_phase.pvar | head -400 > "$WDIR/w4.txt"

bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w4.txt" 1 2 header-read --score-col-nums 3-6 \
  --give-best-pgs --out "$WDIR/gb_degen" > /dev/null 2>&1

EXPECTED_NF=$(printf '%s' "$EXPECTED_HEADER" | awk -F'\t' '{print NF}')
BAD_NF=$(awk -F'\t' -v want="$EXPECTED_NF" 'NF != want {print NR": "NF" fields"}' \
         "$WDIR/gb_degen.best_pgs")
if [ -z "$BAD_NF" ]; then
  echo "PASS  give-best-pgs: every row has $EXPECTED_NF fields (unconverged included)"
else
  echo "FAIL  give-best-pgs: rows with the wrong field count:"
  printf '%s\n' "$BAD_NF"; exit 1
fi

# The degenerate model must be present, last, and NA from BETA to the end.
DEGEN_ROW=$(awk -F'\t' -v id="$C_PGS_ID" '$id=="W4"{print NR}' \
            "$WDIR/gb_degen.best_pgs")
DEGEN_BAD=$(awk -F'\t' -v id="$C_PGS_ID" -v beta="$C_BETA" \
            '$id=="W4"{for(i=beta;i<=NF;i++) if($i!="NA") print i": "$i}' \
            "$WDIR/gb_degen.best_pgs")
if [ "$DEGEN_ROW" = "5" ] && [ -z "$DEGEN_BAD" ]; then
  echo "PASS  give-best-pgs: unconverged row sorts last and is NA from BETA on"
else
  echo "FAIL  give-best-pgs: unconverged row at line '$DEGEN_ROW' (want 5), non-NA: $DEGEN_BAD"
  exit 1
fi

# Case 1: three data rows
ROW_COUNT=$(grep -v '^#' "$BEST" | wc -l | tr -d ' ')
if [ "$ROW_COUNT" -eq 3 ]; then
  echo "PASS  give-best-pgs: 3 data rows"
else
  echo "FAIL  give-best-pgs: expected 3 data rows, got $ROW_COUNT"; exit 1
fi

# Case 2: R2_NULL identical on every row
UNIQ_NULL=$(grep -v '^#' "$BEST" | awk -v c="$C_R2_NULL" '{print $c}' | sort -u | wc -l | tr -d ' ')
if [ "$UNIQ_NULL" -eq 1 ]; then
  echo "PASS  give-best-pgs: R2_NULL identical on all rows"
else
  NULLS=$(grep -v '^#' "$BEST" | awk -v c="$C_R2_NULL" '{print $c}' | tr '\n' ' ')
  echo "FAIL  give-best-pgs: R2_NULL must be identical on every row (same null model), got: $NULLS"; exit 1
fi

# Case 3: rows sorted by R2_INCR descending
ACTUAL_ORDER=$(grep -v '^#' "$BEST" | awk -v c="$C_R2_INCR" '{print $c}')
SORTED_ORDER=$(printf '%s\n' "$ACTUAL_ORDER" | sort -gr)
if [ "$ACTUAL_ORDER" = "$SORTED_ORDER" ]; then
  echo "PASS  give-best-pgs: rows sorted by R2_INCR descending"
else
  echo "FAIL  give-best-pgs: rows not sorted by R2_INCR descending"
  echo "  expected: $(printf '%s ' $SORTED_ORDER)"
  echo "  got:      $(printf '%s ' $ACTUAL_ORDER)"
  exit 1
fi

# ---------------------------------------------------------------------------
# Case 3b: the interval contains its own point estimate, on every row.
#
# This is the invariant the analytic CI stands or falls on, and it is not free:
# squaring an interval on r that straddles zero has minimum 0, not lo_r^2, so
# the naive square puts the lower bound *above* dR2 exactly when the effect is
# near null -- which is what this QT fixture is made of, and what a CI is for.
# Asserted here on QT and again on the CC run, which covers the other branch
# (an interval strictly to the right of zero).
# ---------------------------------------------------------------------------
QT_BAD_CI=$(grep -v '^#' "$BEST" | \
  awk -v id="$C_PGS_ID" -v l="$C_R2_INCR_L95" -v v="$C_R2_INCR" -v u="$C_R2_INCR_U95" '
    $l=="NA" || $v=="NA" || $u=="NA" { print $id": NA in ["$l", "$v", "$u"]"; next }
    !($l+0 <= $v+0 && $v+0 <= $u+0)  { print $id": ["$l", "$v", "$u"]" }')
if [ -z "$QT_BAD_CI" ]; then
  echo "PASS  give-best-pgs: R2_INCR_L95 <= R2_INCR <= R2_INCR_U95 on every row"
else
  echo "FAIL  give-best-pgs: the R2_INCR interval does not contain R2_INCR:"
  printf '%s\n' "$QT_BAD_CI"; exit 1
fi

# Case 3c: the straddle-zero branch, end to end.  Every model here is noise
# against this phenotype, so the interval on r spans zero and the lower bound on
# dR2 must be exactly 0 -- not the small positive lo_r^2 the naive square gives.
# The P check is half the assertion, not decoration: if a future fixture made
# one of these models significant, "exactly 0" would stop being the expected
# answer, and this says so rather than failing mysteriously.
QT_BAD_ZERO=$(grep -v '^#' "$BEST" | \
  awk -v id="$C_PGS_ID" -v p="$C_P" -v l="$C_R2_INCR_L95" \
      '!($p+0 > 0.4 && $l == "0") { print $id": P="$p" R2_INCR_L95="$l }')
if [ -z "$QT_BAD_ZERO" ]; then
  echo "PASS  give-best-pgs: near-null models get R2_INCR_L95 exactly 0"
else
  echo "FAIL  give-best-pgs: near-null model without a zero lower bound:"
  printf '%s\n' "$QT_BAD_ZERO"; exit 1
fi

# Case 4: PHENO column reads QT1 on every row
BAD_PHENO=$(grep -v '^#' "$BEST" | awk -v c="$C_PHENO" '$c != "QT1" {print $c}')
if [ -z "$BAD_PHENO" ]; then
  echo "PASS  give-best-pgs: PHENO column is QT1 on all rows"
else
  echo "FAIL  give-best-pgs: expected PHENO=QT1, got: $BAD_PHENO"; exit 1
fi

# Case 5: categorical phenotype is refused (use SuperPop from the .psam directly)
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno-name SuperPop \
     --covar "$WDIR/covar.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --out "$WDIR/gb_cat" > "$WDIR/cat_out.txt" 2>&1
then
  echo "FAIL  give-best-pgs: categorical phenotype should have been rejected"; exit 1
else
  if grep -qi "categorical" "$WDIR/cat_out.txt"; then
    echo "PASS  give-best-pgs: categorical phenotype refused with 'categorical' message"
  else
    echo "FAIL  give-best-pgs: command failed but output did not mention 'categorical'"
    cat "$WDIR/cat_out.txt"
    exit 1
  fi
fi

# Case 6: --give-best-pgs without --score is refused
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --give-best-pgs --out "$WDIR/gb_ns" > "$WDIR/ns_out.txt" 2>&1
then
  echo "FAIL  give-best-pgs: should have been refused without --score"; exit 1
else
  if grep -qi "score" "$WDIR/ns_out.txt"; then
    echo "PASS  give-best-pgs: refused without --score"
  else
    echo "FAIL  give-best-pgs: command failed but output did not mention '--score'"
    cat "$WDIR/ns_out.txt"
    exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Binary end-to-end: synthesise a case/control phenotype from a score-driven
# latent variable, then check that --give-best-pgs ranks the driver first.
#
# The critical lesson from this project: raw scores are ~1e-3.  An
# unstandardised latent coefficient is swamped by noise and silently ranks
# the wrong model.  Standardise the driving score before use.
# ---------------------------------------------------------------------------

# Step 1: score the three models (reuse the same weight file).
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --out "$WDIR/pgs_for_cc" > /dev/null 2>&1

# Step 2: standardise the W2 model.  Its .sscore column is found by name too --
# which model drives the phenotype is the whole point of the CC section, so it
# is the last thing that should be a bare index.
C_W2_AVG=$(col "$WDIR/pgs_for_cc.sscore" W2_AVG)

# Compute mean and SD from the .sscore file.
read CC_MEAN CC_SD <<EOF
$(awk -v c="$C_W2_AVG" 'NR>1{s+=$c; ss+=$c*$c; n++} END{m=s/n; printf "%.15g %.15g", m, sqrt(ss/n-m*m)}' \
  "$WDIR/pgs_for_cc.sscore")
EOF

# Step 3: dichotomise: z = (score - mean)/sd; logit = 2*z; case if rand < sigmoid(logit).
# This guarantees W2 is the causal driver and should rank first.
awk -v m="$CC_MEAN" -v sd="$CC_SD" -v c="$C_W2_AVG" \
  'BEGIN{OFS="\t"; srand(77); print "#IID","CC_PHENO"}
   NR>1{z=($c-m)/sd; logit=2*z; p=1/(1+exp(-logit));
        print $1,(rand()<p)?2:1}' \
  "$WDIR/pgs_for_cc.sscore" > "$WDIR/pheno_cc.txt"

# Step 4: run give-best-pgs (no --prevalence).
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno_cc.txt" --pheno-name CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --out "$WDIR/gb_cc" > /dev/null 2>&1

CC_BEST="$WDIR/gb_cc.best_pgs"

# Case 7: three data rows
CC_ROW_COUNT=$(grep -v '^#' "$CC_BEST" | wc -l | tr -d ' ')
if [ "$CC_ROW_COUNT" -eq 3 ]; then
  echo "PASS  give-best-pgs CC: 3 data rows"
else
  echo "FAIL  give-best-pgs CC: expected 3 data rows, got $CC_ROW_COUNT"; exit 1
fi

# Case 8: TYPE is CC on all rows
BAD_TYPE=$(grep -v '^#' "$CC_BEST" | awk -v c="$C_TYPE" '$c != "CC" {print $c}')
if [ -z "$BAD_TYPE" ]; then
  echo "PASS  give-best-pgs CC: TYPE is CC on all rows"
else
  echo "FAIL  give-best-pgs CC: expected TYPE=CC, got: $BAD_TYPE"; exit 1
fi

# Case 9: N identical on all rows
UNIQ_N=$(grep -v '^#' "$CC_BEST" | awk -v c="$C_N" '{print $c}' | sort -u | wc -l | tr -d ' ')
if [ "$UNIQ_N" -eq 1 ]; then
  echo "PASS  give-best-pgs CC: N identical on all rows"
else
  echo "FAIL  give-best-pgs CC: N not identical on all rows"; exit 1
fi

# Case 10: N_CASE identical on all rows
UNIQ_NCASE=$(grep -v '^#' "$CC_BEST" | awk -v c="$C_N_CASE" '{print $c}' | sort -u | wc -l | tr -d ' ')
if [ "$UNIQ_NCASE" -eq 1 ]; then
  echo "PASS  give-best-pgs CC: N_CASE identical on all rows"
else
  echo "FAIL  give-best-pgs CC: N_CASE not identical on all rows"; exit 1
fi

# Case 11: R2_INCR sorted descending
CC_ACTUAL_ORDER=$(grep -v '^#' "$CC_BEST" | awk -v c="$C_R2_INCR" '{print $c}')
CC_SORTED_ORDER=$(printf '%s\n' "$CC_ACTUAL_ORDER" | sort -gr)
if [ "$CC_ACTUAL_ORDER" = "$CC_SORTED_ORDER" ]; then
  echo "PASS  give-best-pgs CC: rows sorted by R2_INCR descending"
else
  echo "FAIL  give-best-pgs CC: rows not sorted by R2_INCR descending"
  echo "  expected: $(printf '%s ' $CC_SORTED_ORDER)"
  echo "  got:      $(printf '%s ' $CC_ACTUAL_ORDER)"
  exit 1
fi

# Case 12: the model built from (W2) ranks first.
CC_TOP=$(awk -v c="$C_PGS_ID" 'NR==2{print $c}' "$CC_BEST")
if [ "$CC_TOP" = "W2" ]; then
  echo "PASS  give-best-pgs CC: causal model W2 ranks first"
else
  echo "FAIL  give-best-pgs CC: expected W2 first, got $CC_TOP"; exit 1
fi

# Case 12b: containment again, on a fixture where the interval sits strictly to
# the right of zero (W2 is the causal driver here, R2_INCR ~ 0.35).  Case 3c
# covers the straddle-zero branch; between them both sides of the square are
# exercised end to end.
CC_BAD_CI=$(grep -v '^#' "$CC_BEST" | \
  awk -v id="$C_PGS_ID" -v l="$C_R2_INCR_L95" -v v="$C_R2_INCR" -v u="$C_R2_INCR_U95" '
    $l=="NA" || $v=="NA" || $u=="NA" { print $id": NA in ["$l", "$v", "$u"]"; next }
    !($l+0 <= $v+0 && $v+0 <= $u+0)  { print $id": ["$l", "$v", "$u"]" }')
if [ -z "$CC_BAD_CI" ]; then
  echo "PASS  give-best-pgs CC: R2_INCR_L95 <= R2_INCR <= R2_INCR_U95 on every row"
else
  echo "FAIL  give-best-pgs CC: the R2_INCR interval does not contain R2_INCR:"
  printf '%s\n' "$CC_BAD_CI"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 13: --prevalence actually changes the liability-scale columns.
# Run the same analysis twice, with and without --prevalence 0.05.
# The point estimate and both its bounds are NA together without it and
# numeric together with it -- one column of the three being right proves
# nothing about the other two.
# ---------------------------------------------------------------------------

# Without --prevalence: all three R2L columns are NA on every row.
BAD_R2L_NO=$(grep -v '^#' "$CC_BEST" | \
  awk -v id="$C_PGS_ID" -v v="$C_R2L_INCR" -v l="$C_R2L_INCR_L95" -v u="$C_R2L_INCR_U95" \
      '$v!="NA" || $l!="NA" || $u!="NA" { print $id": "$v" "$l" "$u }')
if [ -z "$BAD_R2L_NO" ]; then
  echo "PASS  give-best-pgs CC: R2L_INCR and both its bounds are NA without --prevalence"
else
  echo "FAIL  give-best-pgs CC: R2L columns should be NA without --prevalence:"
  printf '%s\n' "$BAD_R2L_NO"; exit 1
fi

# With --prevalence 0.05: all three must be numeric on all rows.
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno_cc.txt" --pheno-name CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --prevalence 0.05 --out "$WDIR/gb_cc_prev" > /dev/null 2>&1

BAD_R2L_WITH=$(grep -v '^#' "$WDIR/gb_cc_prev.best_pgs" | \
  awk -v id="$C_PGS_ID" -v v="$C_R2L_INCR" -v l="$C_R2L_INCR_L95" -v u="$C_R2L_INCR_U95" \
      '$v=="NA" || $l=="NA" || $u=="NA" { print $id": "$v" "$l" "$u }')
if [ -z "$BAD_R2L_WITH" ]; then
  echo "PASS  give-best-pgs CC: R2L_INCR and both its bounds are numeric with --prevalence 0.05"
else
  echo "FAIL  give-best-pgs CC: R2L columns should be numeric with --prevalence:"
  printf '%s\n' "$BAD_R2L_WITH"; exit 1
fi

# Case 13b: the liability transform is a multiplication by a positive constant,
# so it must carry containment across unchanged.  If it were ever applied to the
# bounds with the wrong sign or a different prevalence than the point estimate,
# this is where it shows.
BAD_R2L_CI=$(grep -v '^#' "$WDIR/gb_cc_prev.best_pgs" | \
  awk -v id="$C_PGS_ID" -v l="$C_R2L_INCR_L95" -v v="$C_R2L_INCR" -v u="$C_R2L_INCR_U95" \
      '!($l+0 <= $v+0 && $v+0 <= $u+0) { print $id": ["$l", "$v", "$u"]" }')
if [ -z "$BAD_R2L_CI" ]; then
  echo "PASS  give-best-pgs CC: R2L_INCR_L95 <= R2L_INCR <= R2L_INCR_U95 on every row"
else
  echo "FAIL  give-best-pgs CC: the R2L_INCR interval does not contain R2L_INCR:"
  printf '%s\n' "$BAD_R2L_CI"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 14: constant covariate is dropped with a named warning (QT path).
# Build a covariate file with one non-constant column (C1) and one constant
# column (CONST = 5.0 for every sample).
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tC1\tCONST"} NR>1{print $1"\t"(NR%7)/3.0"\t5.0"}' \
  "$PSAM" > "$WDIR/covar_const.txt"

# Run must succeed (exit 0).
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --covar "$WDIR/covar_const.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --out "$WDIR/gb_const" > "$WDIR/const_out.txt" 2>&1
then
  echo "PASS  give-best-pgs: constant covariate does not abort the run"
else
  echo "FAIL  give-best-pgs: constant covariate aborted the run (should have been dropped)"
  cat "$WDIR/const_out.txt"
  exit 1
fi

# Warning must mention 'CONST'.
if grep -qi "CONST" "$WDIR/const_out.txt"; then
  echo "PASS  give-best-pgs: warning names the constant covariate CONST"
else
  echo "FAIL  give-best-pgs: warning did not name the constant covariate"
  cat "$WDIR/const_out.txt"
  exit 1
fi

# Output must have 3 data rows.
CONST_ROW_COUNT=$(grep -v '^#' "$WDIR/gb_const.best_pgs" | wc -l | tr -d ' ')
if [ "$CONST_ROW_COUNT" -eq 3 ]; then
  echo "PASS  give-best-pgs: constant-covariate run produces 3 data rows"
else
  echo "FAIL  give-best-pgs: constant-covariate run produced $CONST_ROW_COUNT rows (expected 3)"; exit 1
fi

# R2_NULL must be identical on every row (null model is still well-defined).
CONST_UNIQ_NULL=$(grep -v '^#' "$WDIR/gb_const.best_pgs" | awk -v c="$C_R2_NULL" '{print $c}' | sort -u | wc -l | tr -d ' ')
if [ "$CONST_UNIQ_NULL" -eq 1 ]; then
  echo "PASS  give-best-pgs: R2_NULL identical on all rows after constant-covariate drop"
else
  CONST_NULLS=$(grep -v '^#' "$WDIR/gb_const.best_pgs" | awk -v c="$C_R2_NULL" '{print $c}' | tr '\n' ' ')
  echo "FAIL  give-best-pgs: R2_NULL not identical after constant-covariate drop: $CONST_NULLS"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 15: all covariates constant — intercept-only fallback, still succeeds.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tA\tB"} NR>1{print $1"\t3.0\t7.0"}' \
  "$PSAM" > "$WDIR/covar_all_const.txt"

if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --covar "$WDIR/covar_all_const.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --out "$WDIR/gb_all_const" > "$WDIR/all_const_out.txt" 2>&1
then
  echo "PASS  give-best-pgs: all-constant covariates falls back to intercept-only"
else
  echo "FAIL  give-best-pgs: all-constant covariates aborted the run"
  cat "$WDIR/all_const_out.txt"
  exit 1
fi

ALL_CONST_ROWS=$(grep -v '^#' "$WDIR/gb_all_const.best_pgs" | wc -l | tr -d ' ')
if [ "$ALL_CONST_ROWS" -eq 3 ]; then
  echo "PASS  give-best-pgs: all-constant intercept-only run produces 3 data rows"
else
  echo "FAIL  give-best-pgs: all-constant intercept-only run produced $ALL_CONST_ROWS rows (expected 3)"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 16-21: multi-phenotype run.
# Build a combined pheno file with one QT and one CC column, use W2 as the
# causal CC driver (same construction as the CC section above), then run
# --pheno-name QT1 CC_PHENO to select both.
# ---------------------------------------------------------------------------

# Reuse the .sscore already produced in the CC section above.
awk -v m="$CC_MEAN" -v sd="$CC_SD" -v c="$C_W2_AVG" \
  'BEGIN{OFS="\t"; srand(77); print "#IID\tQT1\tCC_PHENO"}
   NR>1{iid=$1; qt=(NR%97)/13.0; z=($c-m)/sd; logit=2*z; p=1/(1+exp(-logit));
        cc=(rand()<p)?2:1; print iid, qt, cc}' \
  "$WDIR/pgs_for_cc.sscore" > "$WDIR/multi_pheno.txt"

bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/multi_pheno.txt" --pheno-name QT1 CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --out "$WDIR/gb_multi" > /dev/null 2>&1

MULTI_BEST="$WDIR/gb_multi.best_pgs"

# Case 16: exactly 6 data rows (3 models x 2 phenotypes).
MULTI_ROW_COUNT=$(grep -v '^#' "$MULTI_BEST" | wc -l | tr -d ' ')
if [ "$MULTI_ROW_COUNT" -eq 6 ]; then
  echo "PASS  give-best-pgs multi: 6 data rows (3 models x 2 phenotypes)"
else
  echo "FAIL  give-best-pgs multi: expected 6 data rows, got $MULTI_ROW_COUNT"; exit 1
fi

# Case 17: both phenotype names appear in the PHENO column.
HAS_QT1=$(grep -v '^#' "$MULTI_BEST" | awk -v c="$C_PHENO" '{print $c}' | grep -c '^QT1$')
HAS_CC=$(grep -v '^#' "$MULTI_BEST" | awk -v c="$C_PHENO" '{print $c}' | grep -c '^CC_PHENO$')
if [ "$HAS_QT1" -eq 3 ] && [ "$HAS_CC" -eq 3 ]; then
  echo "PASS  give-best-pgs multi: both PHENO values present (3 rows each)"
else
  echo "FAIL  give-best-pgs multi: expected 3 QT1 and 3 CC_PHENO rows, got QT1=$HAS_QT1 CC_PHENO=$HAS_CC"; exit 1
fi

# Case 18: TYPE is QT for QT1 rows and CC for CC_PHENO rows.
BAD_TYPES=$(grep -v '^#' "$MULTI_BEST" | awk -v p="$C_PHENO" -v t="$C_TYPE" \
  '($p=="QT1" && $t!="QT") || ($p=="CC_PHENO" && $t!="CC") {print $p,$t}')
if [ -z "$BAD_TYPES" ]; then
  echo "PASS  give-best-pgs multi: TYPE is QT for QT1 and CC for CC_PHENO"
else
  echo "FAIL  give-best-pgs multi: unexpected TYPE values: $BAD_TYPES"; exit 1
fi

# Case 19: R2_INCR is descending WITHIN each phenotype group (not globally).
# Check QT1 group (rows 1-3) and CC_PHENO group (rows 4-6) separately.
QT_R2_INCR=$(grep -v '^#' "$MULTI_BEST" | awk -v p="$C_PHENO" -v c="$C_R2_INCR" '$p=="QT1"{print $c}')
CC_R2_INCR=$(grep -v '^#' "$MULTI_BEST" | awk -v p="$C_PHENO" -v c="$C_R2_INCR" '$p=="CC_PHENO"{print $c}')
QT_SORTED=$(printf '%s\n' "$QT_R2_INCR" | sort -gr)
CC_SORTED=$(printf '%s\n' "$CC_R2_INCR" | sort -gr)
if [ "$QT_R2_INCR" = "$QT_SORTED" ] && [ "$CC_R2_INCR" = "$CC_SORTED" ]; then
  echo "PASS  give-best-pgs multi: R2_INCR descending within each phenotype group"
else
  echo "FAIL  give-best-pgs multi: R2_INCR not sorted within phenotype groups"
  echo "  QT R2_INCR:  $(printf '%s ' $QT_R2_INCR)"
  echo "  CC R2_INCR:  $(printf '%s ' $CC_R2_INCR)"
  exit 1
fi

# Case 20: QT rows have NA in OR.
QT_BAD_OR=$(grep -v '^#' "$MULTI_BEST" | awk -v p="$C_PHENO" -v o="$C_OR" '$p=="QT1" && $o!="NA"{print $o}')
if [ -z "$QT_BAD_OR" ]; then
  echo "PASS  give-best-pgs multi: QT rows have NA in OR column"
else
  echo "FAIL  give-best-pgs multi: QT rows should have NA in OR, got: $QT_BAD_OR"; exit 1
fi

# Case 21: QT rows have NA in all three AUC columns.
QT_BAD_AUC=$(grep -v '^#' "$MULTI_BEST" | \
  awk -v p="$C_PHENO" -v an="$C_AUC_NULL" -v af="$C_AUC_FULL" -v ai="$C_AUC_INCR" \
      '$p=="QT1" && ($an!="NA" || $af!="NA" || $ai!="NA"){print $an,$af,$ai}')
if [ -z "$QT_BAD_AUC" ]; then
  echo "PASS  give-best-pgs multi: QT rows have NA in AUC columns"
else
  echo "FAIL  give-best-pgs multi: QT rows should have NA in AUC columns, got: $QT_BAD_AUC"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 21b: the two "_INCR" columns really are the differences their names and
# the documentation claim -- R2_INCR = R2_FULL - R2_NULL on every row, and
# AUC_INCR = AUC_FULL - AUC_NULL on the CC rows that have an AUC at all.
#
# This is the file that has both trait types in it, so one case covers both.
# It is also what stops R2_FULL and the three AUC columns from standing in for
# each other: without it, every assertion about them is "identical on all rows"
# or "is NA", which their neighbours satisfy just as well.  Tolerance is
# relative because these are the printed values, good to six significant
# digits, not the doubles they came from.
# ---------------------------------------------------------------------------
BAD_INCR=$(grep -v '^#' "$MULTI_BEST" | \
  awk -v id="$C_PGS_ID" -v t="$C_TYPE" \
      -v rn="$C_R2_NULL" -v rf="$C_R2_FULL" -v ri="$C_R2_INCR" \
      -v an="$C_AUC_NULL" -v af="$C_AUC_FULL" -v ai="$C_AUC_INCR" '
    function abs(x) { return x < 0 ? -x : x }
    function off(lo, hi, d,    e) {
      e = abs((hi - lo) - d)
      return e > 1e-4 * (abs(lo) + abs(hi) + abs(d))
    }
    off($rn, $rf, $ri) { print $id" R2: "$rf" - "$rn" != "$ri }
    $t=="CC" && off($an, $af, $ai) { print $id" AUC: "$af" - "$an" != "$ai }')
if [ -z "$BAD_INCR" ]; then
  echo "PASS  give-best-pgs multi: R2_INCR and AUC_INCR are the differences they claim"
else
  echo "FAIL  give-best-pgs multi: an _INCR column is not the difference of its pair:"
  printf '%s\n' "$BAD_INCR"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 22-30: --sscore standalone input mode.
# The strongest check: the direct and via-file routes must be identical.
# ---------------------------------------------------------------------------

# Run the normal path to produce a .sscore, then feed it back through --sscore.
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/direct" > /dev/null 2>&1

# Inspect the header to find which columns are scores.
# Our writer emits: #IID QT1 ALLELE_CT NAMED_ALLELE_DOSAGE_SUM W1_AVG W2_AVG W3_AVG
# Score columns are the _AVG columns; determine their 1-based positions.
SCORE_COL_RANGE=$(awk 'NR==1{
  # Strip leading #
  if (substr($0,1,1)=="#") $0=substr($0,2)
  n=split($0,a,"\t")
  first=-1; last=-1
  for(i=1;i<=n;i++){
    if (a[i] ~ /_AVG$|_SUM$/ && a[i] != "NAMED_ALLELE_DOSAGE_SUM") {
      if(first<0) first=i
      last=i
    }
  }
  print first"-"last
  exit
}' "$WDIR/direct.sscore")

# Case 22: --sscore + --score-col-nums produces the same .best_pgs as the direct route.
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/viafile" > /dev/null 2>&1

grep -v '^#' "$WDIR/direct.best_pgs" > "$WDIR/direct_data.txt"
grep -v '^#' "$WDIR/viafile.best_pgs" > "$WDIR/viafile_data.txt"
if diff "$WDIR/direct_data.txt" "$WDIR/viafile_data.txt" > /dev/null 2>&1; then
  echo "PASS  --sscore: direct and via-file routes produce identical output"
else
  echo "FAIL  --sscore: direct and via-file routes differ:"
  diff "$WDIR/direct_data.txt" "$WDIR/viafile_data.txt"
  exit 1
fi

# Case 23: --sscore produces the expected number of data rows.
SSCORE_ROW_COUNT=$(grep -v '^#' "$WDIR/viafile.best_pgs" | wc -l | tr -d ' ')
if [ "$SSCORE_ROW_COUNT" -eq 3 ]; then
  echo "PASS  --sscore: 3 data rows"
else
  echo "FAIL  --sscore: expected 3 data rows, got $SSCORE_ROW_COUNT"; exit 1
fi

# Case 24: subset test — a .sscore whose sample set is a strict subset of the
# phenotype file.  Drop the last 100 samples from the .sscore; warn+run.
HEAD_LINES=$(($(wc -l < "$WDIR/direct.sscore") - 100))
awk "NR<=1 || NR<=$HEAD_LINES" "$WDIR/direct.sscore" > "$WDIR/subset.sscore"

SUBSET_OUT=$(bin/adjutant \
  --sscore "$WDIR/subset.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/subset_result" 2>&1)

if echo "$SUBSET_OUT" | grep -qi "dropped"; then
  echo "PASS  --sscore: subset .sscore emits a drop warning"
else
  echo "FAIL  --sscore: subset .sscore should have emitted a drop warning"
  echo "$SUBSET_OUT"
  exit 1
fi

# The N in the subset output should be < 2504.
SUBSET_N=$(grep -v '^#' "$WDIR/subset_result.best_pgs" | awk -v c="$C_N" '{print $c}' | sort -u)
if [ "$SUBSET_N" -lt 2504 ]; then
  echo "PASS  --sscore: N in subset output reflects the intersection ($SUBSET_N < 2504)"
else
  echo "FAIL  --sscore: expected N < 2504 in subset output, got $SUBSET_N"; exit 1
fi

# Case 25: --sscore combined with --pfile is refused.
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --sscore "$WDIR/direct.sscore" \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --give-best-pgs --out "$WDIR/conflict_test" > "$WDIR/conflict_out.txt" 2>&1
then
  echo "FAIL  --sscore: should have been refused when combined with --pfile"; exit 1
else
  if grep -qi "pfile\|conflict\|cannot be used\|mutually exclusive" "$WDIR/conflict_out.txt"; then
    echo "PASS  --sscore: refused when combined with --pfile"
  else
    echo "FAIL  --sscore: refused but wrong error message"
    cat "$WDIR/conflict_out.txt"
    exit 1
  fi
fi

# Case 26: --sscore without --pheno is refused.
if bin/adjutant \
     --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
     --give-best-pgs --out "$WDIR/nopheno_test" > "$WDIR/nopheno_out.txt" 2>&1
then
  echo "FAIL  --sscore: should have been refused without --pheno"; exit 1
else
  if grep -qi "pheno" "$WDIR/nopheno_out.txt"; then
    echo "PASS  --sscore: refused without --pheno"
  else
    echo "FAIL  --sscore: refused but wrong error message"
    cat "$WDIR/nopheno_out.txt"
    exit 1
  fi
fi

# Case 27: duplicate IIDs in the .sscore are refused.
# Create a .sscore with a duplicate IID (copy the first data row with same IID).
{ head -2 "$WDIR/direct.sscore"; head -2 "$WDIR/direct.sscore" | tail -1; tail -n +3 "$WDIR/direct.sscore"; } \
  > "$WDIR/dupid.sscore"

if bin/adjutant \
     --sscore "$WDIR/dupid.sscore" --score-col-nums "$SCORE_COL_RANGE" \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --give-best-pgs --out "$WDIR/dup_test" > "$WDIR/dup_out.txt" 2>&1
then
  echo "FAIL  --sscore: should have been refused with duplicate IIDs"; exit 1
else
  if grep -qi "duplicate" "$WDIR/dup_out.txt"; then
    echo "PASS  --sscore: duplicate IIDs in .sscore are refused"
  else
    echo "FAIL  --sscore: refused but wrong error message"
    cat "$WDIR/dup_out.txt"
    exit 1
  fi
fi

# Case 28: CC phenotype via --sscore matches the direct route.
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --pheno "$WDIR/pheno_cc.txt" --pheno-name CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/direct_cc" > /dev/null 2>&1

CC_SSCORE_FILE="$WDIR/direct_cc.sscore"
CC_SCORE_RANGE=$(awk 'NR==1{
  if (substr($0,1,1)=="#") $0=substr($0,2)
  n=split($0,a,"\t"); first=-1; last=-1
  for(i=1;i<=n;i++){
    if (a[i] ~ /_AVG$|_SUM$/ && a[i] != "NAMED_ALLELE_DOSAGE_SUM"){
      if(first<0) first=i; last=i
    }
  }
  print first"-"last; exit
}' "$CC_SSCORE_FILE")

bin/adjutant \
  --sscore "$CC_SSCORE_FILE" --score-col-nums "$CC_SCORE_RANGE" \
  --pheno "$WDIR/pheno_cc.txt" --pheno-name CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/viafile_cc" > /dev/null 2>&1

grep -v '^#' "$WDIR/direct_cc.best_pgs" > "$WDIR/direct_cc_data.txt"
grep -v '^#' "$WDIR/viafile_cc.best_pgs" > "$WDIR/viafile_cc_data.txt"
if diff "$WDIR/direct_cc_data.txt" "$WDIR/viafile_cc_data.txt" > /dev/null 2>&1; then
  echo "PASS  --sscore CC: direct and via-file routes produce identical output"
else
  echo "FAIL  --sscore CC: direct and via-file routes differ:"
  diff "$WDIR/direct_cc_data.txt" "$WDIR/viafile_cc_data.txt"
  exit 1
fi

# ---------------------------------------------------------------------------
# Cases 29-41: --bim / --pgs-model end-to-end.
#
# Source file is derived from Tests/TEST_PHASED_VCF/mt.bim with deliberate
# perturbations so every status is represented:
#   rows 1-5   exact      (use bim's A1 A2 directly)
#   rows 6-8   reversed   (swap A1/A2 columns)
#   rows 9-11  flipped    (complement both alleles; only single-char variants)
#   rows 12-14 flipped+reversed (complement and swap)
#   ambiguous  left as-is (the 17 A/T or C/G variants in the bim)
#   rows in    3 positions not in the bim (novel positions)
#   1 position present but alleles irreconcilable
#
# The bim has 109 variants: 17 ambiguous + 92 non-ambiguous single-char pairs
# plus a few indels.
#
# We construct the source file from awk so no Python/R is needed.
# ---------------------------------------------------------------------------
# NOTE: never redirect a run's stdout onto "$WDIR/<out>.log".  adjutant writes
# that file itself from --out, so the shell capture and adjutant's logger become
# two writers on one path at independent offsets.  The result looks fine until
# the byte offsets shift, at which point the file gains a second copy of the
# summary and any "grep -Eo" extraction silently returns two numbers.  Capture
# stdout as <out>_cli.txt instead.
BIM="Tests/TEST_PHASED_VCF/mt.bim"
PFX="Tests/TEST_PHASED_VCF/mt_mixed_phase"

# Build source file.
# Strategy: iterate over the bim in awk.
# - Track which variants are single-char non-ambiguous (can flip).
# - Assign a perturbation: exact, reversed, flipped, flipped+rev, or ambiguous.
# - Also inject 3 phantom positions and 1 bad-allele row.
#
# Complement function in awk: comp[A]=T etc.
# .bim columns: $1=CHR $2=ID $3=CM $4=POS $5=A1 $6=A2
# Source file columns: CHR POS A1(effect) A2(other) BETA
# In source convention: A1=effect=a1, A2=other=a0.
# In bim convention:    col5=A1=a1(effect), col6=A2=a0(other).
# So for EXACT match: src_a1=bim_A1, src_a0=bim_A2.
# For REVERSED: src_a1=bim_A2, src_a0=bim_A1.
awk '
BEGIN {
  OFS="\t"
  comp["A"]="T"; comp["T"]="A"; comp["C"]="G"; comp["G"]="C"
  exact_ct=0; rev_ct=0; flip_ct=0; fliprev_ct=0; amb_ct=0
  inj=0
}
{
  chr=$1; pos=$4; a1=$5; a0=$6
  # Determine if single-char and ambiguous
  single = (length(a1)==1 && length(a0)==1)
  ambig  = single && ((a1=="A"&&a0=="T")||(a1=="T"&&a0=="A")||(a1=="C"&&a0=="G")||(a1=="G"&&a0=="C"))

  if (ambig) {
    # Leave ambiguous as-is: src_a1=a1, src_a0=a0
    print chr, pos, a1, a0, 0.05
    amb_ct++
  } else if (exact_ct < 5) {
    # EXACT: src_a1=a1, src_a0=a0
    print chr, pos, a1, a0, 0.10
    exact_ct++
  } else if (rev_ct < 3) {
    # REVERSED: swap
    print chr, pos, a0, a1, 0.20
    rev_ct++
  } else if (flip_ct < 3 && single && comp[a1] != "" && comp[a0] != "") {
    # FLIPPED: complement both
    print chr, pos, comp[a1], comp[a0], 0.30
    flip_ct++
  } else if (fliprev_ct < 3 && single && comp[a1] != "" && comp[a0] != "") {
    # FLIPPED+REVERSED: complement and swap
    print chr, pos, comp[a0], comp[a1], 0.40
    fliprev_ct++
  } else {
    # Remaining: use exact
    print chr, pos, a1, a0, 0.05
  }
}
END {
  # 3 phantom positions (not in bim)
  print "1", "9999999", "A", "G", 0.50
  print "2", "9999998", "C", "T", 0.50
  print "3", "9999997", "G", "A", 0.50
  # 1 bad allele row (pos in bim but alleles irreconcilable)
  # 1_9894436_C_T has A1=T A2=C; use X/Y which cannot match
  print "1", "9894436", "X", "Y", 0.99
}
' "$BIM" > "$WDIR/pgs_src.txt"

# Case 29: run --pgs-model without --allow-ambiguous.
bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src.txt" 1 2 3 4 5 \
  --out "$WDIR/harm" > "$WDIR/harm_cli.txt" 2>&1

# The .match_report columns are resolved by name for the same reason the
# .best_pgs ones are: STATUS and OUT_BETA are adjacent to columns of the same
# shape, and a literal index cannot tell them apart.
C_STATUS=$(col "$WDIR/harm.match_report" STATUS)
C_OUT_BETA=$(col "$WDIR/harm.match_report" OUT_BETA)
C_SRC_EFFECT=$(col "$WDIR/harm.match_report" SRC_EFFECT)
C_OUT_A1=$(col "$WDIR/harm.match_report" OUT_A1)

# Number of lines in match_report (excluding header) must equal number of
# source rows (including the 4 injected rows).
SRC_ROWS=$(grep -c '' "$WDIR/pgs_src.txt")
REPORT_ROWS=$(grep -v '^#' "$WDIR/harm.match_report" | wc -l | tr -d ' ')
if [ "$REPORT_ROWS" -eq "$SRC_ROWS" ]; then
  echo "PASS  --pgs-model: match_report has one row per source variant"
else
  echo "FAIL  --pgs-model: expected $SRC_ROWS rows in report, got $REPORT_ROWS"
  exit 1
fi

# Case 30: REVERSED variants must have NEGATED beta in .matched.
# We made reversed variants (rows 6-8) with src_beta=0.20.
# Their out_beta must be -0.20.
BAD_REV=$(awk -v st="$C_STATUS" -v b="$C_OUT_BETA" \
           'NR>1 && $st=="REVERSED" && $b > -0.19 {print NR, $0}' \
           "$WDIR/harm.match_report" | head -3)
if [ -z "$BAD_REV" ]; then
  echo "PASS  --pgs-model: REVERSED variants have negated beta"
else
  echo "FAIL  --pgs-model: REVERSED variant with wrong beta:"
  echo "$BAD_REV"
  exit 1
fi

# Case 31: EXACT variants must have unchanged beta.
BAD_EXACT=$(awk -v st="$C_STATUS" -v b="$C_OUT_BETA" \
  'NR>1 && $st=="EXACT" && ($b != "0.1" && $b != "0.05") {
  # Exact cases have beta 0.10 or 0.05
  if ($b < 0) print NR, $0
}' "$WDIR/harm.match_report" | head -3)
if [ -z "$BAD_EXACT" ]; then
  echo "PASS  --pgs-model: EXACT variants have unchanged (positive) beta"
else
  echo "FAIL  --pgs-model: EXACT variant with wrong beta: $BAD_EXACT"
  exit 1
fi

# Case 32: 3 phantom positions must appear as NO_MATCH_POS in report.
NO_POS_CT=$(awk -v st="$C_STATUS" '$st=="NO_MATCH_POS"' "$WDIR/harm.match_report" | wc -l | tr -d ' ')
if [ "$NO_POS_CT" -ge 3 ]; then
  echo "PASS  --pgs-model: phantom positions appear as NO_MATCH_POS ($NO_POS_CT rows)"
else
  echo "FAIL  --pgs-model: expected >=3 NO_MATCH_POS rows, got $NO_POS_CT"
  exit 1
fi

# Case 33: 1 bad-allele row must appear as NO_MATCH_ALLELE or DUPLICATE.
NO_ALLELE_CT=$(awk -v st="$C_STATUS" '$st=="NO_MATCH_ALLELE" || $st=="DUPLICATE"' \
               "$WDIR/harm.match_report" | wc -l | tr -d ' ')
if [ "$NO_ALLELE_CT" -ge 1 ]; then
  echo "PASS  --pgs-model: irreconcilable allele appears as NO_MATCH_ALLELE/DUPLICATE"
else
  echo "FAIL  --pgs-model: expected >=1 NO_MATCH_ALLELE/DUPLICATE row, got $NO_ALLELE_CT"
  exit 1
fi

# Case 34: ambiguous variants (17) must appear as AMBIGUOUS in report.
AMB_CT=$(awk -v st="$C_STATUS" '$st=="AMBIGUOUS"' "$WDIR/harm.match_report" | wc -l | tr -d ' ')
if [ "$AMB_CT" -eq 17 ]; then
  echo "PASS  --pgs-model: 17 ambiguous variants in report"
else
  echo "FAIL  --pgs-model: expected 17 AMBIGUOUS rows, got $AMB_CT"
  exit 1
fi

# Case 35: summary line mentions ambiguous count.
if saidflat "$WDIR/harm_cli.txt" "17 ambiguous"; then
  echo "PASS  --pgs-model: summary reports 17 ambiguous"
else
  echo "FAIL  --pgs-model: summary does not mention '17 ambiguous'"
  grep "pgs-model:" "$WDIR/harm_cli.txt"
  exit 1
fi

# Case 36: --allow-ambiguous increases matched count by exactly 17.
MATCHED_NO_AMB=$(grep "matched" "$WDIR/harm_cli.txt" | grep -Eo '[0-9]+ matched' | awk '{print $1}')

bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src.txt" 1 2 3 4 5 \
  --allow-ambiguous \
  --out "$WDIR/harm_amb" > "$WDIR/harm_amb_cli.txt" 2>&1

MATCHED_WITH_AMB=$(grep "matched" "$WDIR/harm_amb_cli.txt" | grep -Eo '[0-9]+ matched' | awk '{print $1}')
DIFF_AMB=$((MATCHED_WITH_AMB - MATCHED_NO_AMB))
if [ "$DIFF_AMB" -eq 17 ]; then
  echo "PASS  --pgs-model: --allow-ambiguous increases matched count by 17"
else
  echo "FAIL  --pgs-model: expected +17 from --allow-ambiguous, got +$DIFF_AMB (base $MATCHED_NO_AMB, with $MATCHED_WITH_AMB)"
  exit 1
fi

# Case 37: .matched file has a header row and only matched variants.
MATCHED_HEADER=$(head -1 "$WDIR/harm.matched")
if [ "$MATCHED_HEADER" = "$(printf 'ID\tA1\tBETA')" ]; then
  echo "PASS  --pgs-model: .matched has correct header"
else
  echo "FAIL  --pgs-model: .matched header is '$MATCHED_HEADER'"
  exit 1
fi

# Case 38: FLIPPED variants have unchanged (positive) beta in .matched.
BAD_FLIP=$(awk -v st="$C_STATUS" -v b="$C_OUT_BETA" \
           'NR>1 && $st=="FLIPPED" && $b < 0 {print NR, $0}' \
           "$WDIR/harm.match_report" | head -3)
if [ -z "$BAD_FLIP" ]; then
  echo "PASS  --pgs-model: FLIPPED variants have unchanged (positive) beta"
else
  echo "FAIL  --pgs-model: FLIPPED variant with negated beta: $BAD_FLIP"
  exit 1
fi

# Case 39: FLIPPED_REVERSED variants have negated beta.
BAD_FR=$(awk -v st="$C_STATUS" -v b="$C_OUT_BETA" \
         'NR>1 && $st=="FLIPPED_REVERSED" && $b >= 0 {print NR, $0}' \
         "$WDIR/harm.match_report" | head -3)
if [ -z "$BAD_FR" ]; then
  echo "PASS  --pgs-model: FLIPPED_REVERSED variants have negated beta"
else
  echo "FAIL  --pgs-model: FLIPPED_REVERSED variant with wrong beta: $BAD_FR"
  exit 1
fi

# Case 40: .matched feeds into --score and produces a .sscore.
bin/adjutant \
  --pfile "$PFX" vzs \
  --score "$WDIR/harm.matched" 1 2 header-read ignore-dup-ids \
  --out "$WDIR/scored_from_harm" > "$WDIR/scored_from_harm_cli.txt" 2>&1

if [ -s "$WDIR/scored_from_harm.sscore" ]; then
  echo "PASS  --pgs-model: .matched feeds directly into --score producing .sscore"
else
  echo "FAIL  --pgs-model: .matched did not produce a .sscore"
  cat "$WDIR/scored_from_harm_cli.txt"
  exit 1
fi

# Case 41: --pgs-model refuses --score combination.
if bin/adjutant \
     --bim "$BIM" \
     --pgs-model "$WDIR/pgs_src.txt" 1 2 3 4 5 \
     --score "$WDIR/pgs_src.txt" 1 2 3 \
     --out "$WDIR/conflict_pgs" > "$WDIR/conflict_pgs_cli.txt" 2>&1
then
  echo "FAIL  --pgs-model: should have been refused with --score"; exit 1
else
  if grep -qi "cannot be used\|requires\|mutually exclusive" "$WDIR/conflict_pgs_cli.txt"; then
    echo "PASS  --pgs-model: refused when combined with --score"
  else
    echo "FAIL  --pgs-model: refused but wrong error message"
    cat "$WDIR/conflict_pgs_cli.txt"
    exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Cases 41a-41m: --resolve-ambiguous, palindromes recovered from allele frequencies.
#
# A second source file, carrying a sixth column (the effect-allele frequency in
# the source study).  The target frequencies are synthesised from the same .bim
# so every number below is known exactly rather than looked up:
#
#   tgt.afreq       REF = .bim A2, ALT = .bim A1, ALT_FREQS = 0.85
#                   => the frequency of the .bim's A1 is 0.85 for every variant
#   tgt_flip.afreq  REF = .bim A1, ALT = .bim A2, ALT_FREQS = 0.15
#                   => the SAME frequencies, written the other way round
#
# The 17 palindromic .bim variants cycle through three source shapes:
#   i%3==0  effect = .bim A1, source freq 0.85  -> RESOLVED_EXACT,        +0.30
#   i%3==1  effect = .bim A2, source freq 0.15  -> RESOLVED_REVERSED,     -0.30
#   i%3==2  effect = .bim A1, source freq 0.48  -> AMBIGUOUS_UNRESOLVED
# giving 6 + 6 + 5.  The i%3==1 rows are the ones the feature exists for: the
# source named the other allele of the pair, which no allele comparison can
# detect, and a backwards frequency comparison would emit +0.30.
#
# One extra row is an A/T source against a C/G target: palindromic on both
# sides, no shared allele, must come out NO_MATCH_ALLELE rather than resolved.
# ---------------------------------------------------------------------------
awk -v OFS='\t' '
{
  chr=$1; pos=$4; a1=$5; a0=$6
  single = (length(a1)==1 && length(a0)==1)
  pal = single && ((a1=="A"&&a0=="T")||(a1=="T"&&a0=="A")|| \
                   (a1=="C"&&a0=="G")||(a1=="G"&&a0=="C"))
  if (!pal) { print chr, pos, a1, a0, 0.10, 0.85; next }
  m = pal_i % 3; pal_i++
  if      (m == 0) print chr, pos, a1, a0, 0.30, 0.85
  else if (m == 1) print chr, pos, a0, a1, 0.30, 0.15
  else             print chr, pos, a1, a0, 0.30, 0.48
}
END { print "2", "182791973", "A", "T", 0.30, 0.85 }
' "$BIM" > "$WDIR/pgs_src_freq.txt"

awk -v OFS='\t' 'BEGIN{print "#CHROM","ID","REF","ALT","ALT_FREQS","OBS_CT"}
                 {print $1,$2,$6,$5,0.85,1000}' "$BIM" > "$WDIR/tgt.afreq"
awk -v OFS='\t' 'BEGIN{print "#CHROM","ID","REF","ALT","ALT_FREQS","OBS_CT"}
                 {print $1,$2,$5,$6,0.15,1000}' "$BIM" > "$WDIR/tgt_flip.afreq"

# Case 41a: the primary regression guard.  Same file, sixth column present,
# --resolve-ambiguous absent: every palindrome is still dropped as AMBIGUOUS and
# nothing is resolved.  Default behaviour must not move.
bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 6 \
  --out "$WDIR/ra_off" > "$WDIR/ra_off_cli.txt" 2>&1

RA_OFF_AMB=$(awk -v st="$C_STATUS" '$st=="AMBIGUOUS"' "$WDIR/ra_off.match_report" | wc -l | tr -d ' ')
RA_OFF_RES=$(awk -v st="$C_STATUS" '$st ~ /^RESOLVED_/ || $st=="AMBIGUOUS_UNRESOLVED"' \
             "$WDIR/ra_off.match_report" | wc -l | tr -d ' ')
if [ "$RA_OFF_AMB" -eq 18 ] && [ "$RA_OFF_RES" -eq 0 ]; then
  echo "PASS  --resolve-ambiguous: absent, all 18 palindromic rows still AMBIGUOUS"
else
  echo "FAIL  --resolve-ambiguous: absent, expected 18 AMBIGUOUS and 0 resolved, got $RA_OFF_AMB and $RA_OFF_RES"
  exit 1
fi

# Case 41b: with --resolve-ambiguous the matched count rises by exactly 12.
bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/tgt.afreq" --resolve-ambiguous \
  --out "$WDIR/ra_on" > "$WDIR/ra_on_cli.txt" 2>&1

RA_OFF_MATCHED=$(grep -Eo '[0-9]+ matched' "$WDIR/ra_off_cli.txt" | head -1 | awk '{print $1}')
RA_ON_MATCHED=$(grep -Eo '[0-9]+ matched' "$WDIR/ra_on_cli.txt" | head -1 | awk '{print $1}')
if [ "$((RA_ON_MATCHED - RA_OFF_MATCHED))" -eq 12 ]; then
  echo "PASS  --resolve-ambiguous: recovers 12 palindromes ($RA_OFF_MATCHED -> $RA_ON_MATCHED matched)"
else
  echo "FAIL  --resolve-ambiguous: expected +12 matched, got $RA_OFF_MATCHED -> $RA_ON_MATCHED"
  exit 1
fi

# Case 41c: 6 of each resolved status, which is how the source file was built.
RE_CT=$(awk -v st="$C_STATUS" '$st=="RESOLVED_EXACT"' "$WDIR/ra_on.match_report" | wc -l | tr -d ' ')
RR_CT=$(awk -v st="$C_STATUS" '$st=="RESOLVED_REVERSED"' "$WDIR/ra_on.match_report" | wc -l | tr -d ' ')
if [ "$RE_CT" -eq 6 ] && [ "$RR_CT" -eq 6 ]; then
  echo "PASS  --resolve-ambiguous: 6 RESOLVED_EXACT and 6 RESOLVED_REVERSED"
else
  echo "FAIL  --resolve-ambiguous: expected 6 and 6 resolved, got $RE_CT and $RR_CT"
  exit 1
fi

# Case 41d: the load-bearing one -- did the frequency comparison recover the
# right orientation, on every row?
#
# OUT_A1 is the .bim's A1.  The resolution never looks at the allele labels
# (for a palindrome it cannot: both hypotheses produce the same labels), so
# whether it lands on the orientation the labels encode is a real question.
# Each source row was written on the target's strand and given the frequency
# OF ITS OWN EFFECT ALLELE, so the answer must be yes every time:
#   SRC_EFFECT == OUT_A1  -> RESOLVED_EXACT,    beta kept at  +0.3
#   SRC_EFFECT != OUT_A1  -> RESOLVED_REVERSED, beta negated to -0.3
#
# Verified by mutation: reversing the frequency comparison relabels the rows
# but leaves "RESOLVED_EXACT rows have +0.3" true, so that weaker assertion
# passed against a build that inverted 6 of the 12 effects.  This one fails.
BAD_ORIENT=$(awk -v st="$C_STATUS" -v b="$C_OUT_BETA" -v e="$C_SRC_EFFECT" -v o="$C_OUT_A1" '
  NR>1 && $st ~ /^RESOLVED_/ {
    want_st = ($e == $o) ? "RESOLVED_EXACT" : "RESOLVED_REVERSED"
    want_b  = ($e == $o) ? "0.3" : "-0.3"
    if ($st != want_st || $b != want_b) print NR": "$0"  (want "want_st" "want_b")"
  }' "$WDIR/ra_on.match_report" | head -3)
if [ -z "$BAD_ORIENT" ]; then
  echo "PASS  --resolve-ambiguous: every resolved row recovers its true orientation and sign"
else
  echo "FAIL  --resolve-ambiguous: resolved row with the wrong orientation or sign:"
  echo "$BAD_ORIENT"
  exit 1
fi

# The same negation must reach .matched, which is the file --score consumes.
MATCHED_NEG=$(awk 'NR>1 && $3=="-0.3"' "$WDIR/ra_on.matched" | wc -l | tr -d ' ')
if [ "$MATCHED_NEG" -eq 6 ]; then
  echo "PASS  --resolve-ambiguous: the 6 negated betas reach .matched"
else
  echo "FAIL  --resolve-ambiguous: expected 6 rows with beta -0.3 in .matched, got $MATCHED_NEG"
  exit 1
fi

# Case 41e: source frequency 0.48 stays unresolved whatever the target says.
UNRES_CT=$(awk -v st="$C_STATUS" '$st=="AMBIGUOUS_UNRESOLVED"' "$WDIR/ra_on.match_report" | wc -l | tr -d ' ')
UNRES_IN_MATCHED=$(awk -v st="$C_STATUS" '$st=="AMBIGUOUS_UNRESOLVED" && $NF != "NA"' \
                   "$WDIR/ra_on.match_report" | wc -l | tr -d ' ')
if [ "$UNRES_CT" -eq 5 ] && [ "$UNRES_IN_MATCHED" -eq 0 ]; then
  echo "PASS  --resolve-ambiguous: 5 near-0.5 palindromes AMBIGUOUS_UNRESOLVED and dropped"
else
  echo "FAIL  --resolve-ambiguous: expected 5 unresolved and 0 emitted, got $UNRES_CT and $UNRES_IN_MATCHED"
  exit 1
fi

# Case 41f: an A/T source against a C/G target shares no allele.  Frequencies
# must not be allowed to "resolve" it.
CROSS=$(awk -v st="$C_STATUS" '$1=="2" && $2=="182791973" && $3=="A" {print $st}' \
        "$WDIR/ra_on.match_report")
if [ "$CROSS" = "NO_MATCH_ALLELE" ]; then
  echo "PASS  --resolve-ambiguous: A/T source vs C/G target is NO_MATCH_ALLELE"
else
  echo "FAIL  --resolve-ambiguous: A/T vs C/G came out '$CROSS', want NO_MATCH_ALLELE"
  exit 1
fi

# Case 41g: nothing guarantees the .afreq's REF/ALT orientation matches the
# .bim's A1/A2.  The same frequencies written the other way round must give
# the same answer -- the direct test that orientation is resolved, not assumed.
bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/tgt_flip.afreq" --resolve-ambiguous \
  --out "$WDIR/ra_flip" > "$WDIR/ra_flip_cli.txt" 2>&1
if diff -q "$WDIR/ra_on.match_report" "$WDIR/ra_flip.match_report" > /dev/null; then
  echo "PASS  --resolve-ambiguous: flipped .afreq REF/ALT orientation gives identical output"
else
  echo "FAIL  --resolve-ambiguous: flipped .afreq orientation changed the result"
  diff "$WDIR/ra_on.match_report" "$WDIR/ra_flip.match_report" | head -6
  exit 1
fi

# Case 41h: the summary accounts for every palindrome it decided.
if saidflat "$WDIR/ra_on_cli.txt" "17 ambiguous; 12 resolved by frequency" &&
   saidflat "$WDIR/ra_on_cli.txt" "5 unresolved and dropped"; then
  echo "PASS  --resolve-ambiguous: summary reports 17 ambiguous, 12 resolved, 5 dropped"
else
  echo "FAIL  --resolve-ambiguous: summary line wrong"
  grep "pgs-model:" "$WDIR/ra_on_cli.txt"
  exit 1
fi

# Case 41i: a multi-allelic .afreq row is unavailable, not guessed at.  Rewrite
# the row for one variant that otherwise resolves and watch it stop resolving.
awk -v OFS='\t' '$2=="1_56520159_T_A" {$4=$4",G"; $5=$5",0.01"} {print}' \
  "$WDIR/tgt.afreq" > "$WDIR/tgt_multi.afreq"
bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/tgt_multi.afreq" --resolve-ambiguous \
  --out "$WDIR/ra_multi" > "$WDIR/ra_multi_cli.txt" 2>&1
MULTI_STATUS=$(awk -v st="$C_STATUS" '$1=="1" && $2=="56520159" {print $st}' \
               "$WDIR/ra_multi.match_report")
if [ "$MULTI_STATUS" = "AMBIGUOUS_UNRESOLVED" ]; then
  echo "PASS  --resolve-ambiguous: multi-allelic .afreq row is unavailable, not guessed"
else
  echo "FAIL  --resolve-ambiguous: multi-allelic row gave '$MULTI_STATUS', want AMBIGUOUS_UNRESOLVED"
  exit 1
fi

# Case 41j: a real plink2 .afreq -- multi-allelic rows, repeated IDs and all --
# joined on variant ID.  Source frequencies are read out of it, so every one of
# the .bim's 17 palindromes has a usable pair and all 17 resolve.
awk -v OFS='\t' '
NR==FNR { if (FNR>1 && !($2 in seen)) { seen[$2]=1; ref[$2]=$3; alt[$2]=$4; f[$2]=$5 } ; next }
{ a1=$5; a0=$6
  if (length(a1)!=1 || length(a0)!=1) next
  if (!((a1=="A"&&a0=="T")||(a1=="T"&&a0=="A")|| \
        (a1=="C"&&a0=="G")||(a1=="G"&&a0=="C"))) next
  id=$2
  if (!(id in seen)) next
  if (alt[id]==a1) sf=f[id]; else if (ref[id]==a1) sf=1-f[id]; else next
  print $1, $4, a1, a0, 0.30, sf
}' Tests/score_fixtures/freqs.afreq "$BIM" > "$WDIR/pgs_src_real.txt"
bin/adjutant \
  --bim "$BIM" \
  --pgs-model "$WDIR/pgs_src_real.txt" 1 2 3 4 5 6 \
  --read-freq Tests/score_fixtures/freqs.afreq --resolve-ambiguous \
  --out "$WDIR/ra_real" > "$WDIR/ra_real_cli.txt" 2>&1
REAL_RES=$(awk -v st="$C_STATUS" '$st ~ /^RESOLVED_/' "$WDIR/ra_real.match_report" | wc -l | tr -d ' ')
if [ "$REAL_RES" -eq 17 ]; then
  echo "PASS  --resolve-ambiguous: all 17 of the bundled .bim's palindromes resolve on a real .afreq"
else
  echo "FAIL  --resolve-ambiguous: expected 17 resolved against the real .afreq, got $REAL_RES"
  grep "pgs-model:" "$WDIR/ra_real_cli.txt"
  exit 1
fi

# Case 41k: refused without --read-freq.
if bin/adjutant --bim "$BIM" \
     --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 6 --resolve-ambiguous \
     --out "$WDIR/ra_e1" > "$WDIR/ra_e1_cli.txt" 2>&1
then
  echo "FAIL  --resolve-ambiguous: should have been refused without --read-freq"; exit 1
elif saidflat "$WDIR/ra_e1_cli.txt" "requires --read-freq"; then
  echo "PASS  --resolve-ambiguous: refused without --read-freq"
else
  echo "FAIL  --resolve-ambiguous: refused without --read-freq but wrong message"
  cat "$WDIR/ra_e1_cli.txt"; exit 1
fi

# Case 41l: refused without the sixth column index.
if bin/adjutant --bim "$BIM" \
     --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 \
     --read-freq "$WDIR/tgt.afreq" --resolve-ambiguous \
     --out "$WDIR/ra_e2" > "$WDIR/ra_e2_cli.txt" 2>&1
then
  echo "FAIL  --resolve-ambiguous: should have been refused without the sixth column"; exit 1
elif saidflat "$WDIR/ra_e2_cli.txt" "sixth --pgs-model column index"; then
  echo "PASS  --resolve-ambiguous: refused without the sixth column index"
else
  echo "FAIL  --resolve-ambiguous: refused without the sixth column but wrong message"
  cat "$WDIR/ra_e2_cli.txt"; exit 1
fi

# Case 41m: refused with --allow-ambiguous.  The two express contradictory
# intent, and mixing blindly-kept and frequency-resolved palindromes in one
# output would hide a real difference in evidence.
if bin/adjutant --bim "$BIM" \
     --pgs-model "$WDIR/pgs_src_freq.txt" 1 2 3 4 5 6 \
     --read-freq "$WDIR/tgt.afreq" --resolve-ambiguous --allow-ambiguous \
     --out "$WDIR/ra_e3" > "$WDIR/ra_e3_cli.txt" 2>&1
then
  echo "FAIL  --resolve-ambiguous: should have been refused with --allow-ambiguous"; exit 1
elif saidflat "$WDIR/ra_e3_cli.txt" "cannot be used with --allow-ambiguous"; then
  echo "PASS  --resolve-ambiguous: refused when combined with --allow-ambiguous"
else
  echo "FAIL  --resolve-ambiguous: refused with --allow-ambiguous but wrong message"
  cat "$WDIR/ra_e3_cli.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 41n-41s: the fixture the frequency comparison cannot fake.
#
# Case 41d above checks that every resolved row agrees with its allele labels.
# That fixture writes every source row on the TARGET's strand, so
# "SRC_EFFECT == OUT_A1  <=>  RESOLVED_EXACT" holds by construction: it catches
# a comparison inversion, but an implementation that ignored frequencies
# altogether and compared allele labels would pass it too.  That is the one
# feature this flag adds, so it gets a fixture that cannot be passed that way.
#
# Every target below has A1 at frequency 0.20.  The d_at_* rows are the point:
#
#   d_at_exact_lo  source effect A, other T, source freq 0.20
#   d_at_rev_lo    source effect A, other T, source freq 0.80
#
# Identical allele labels, differing only in the source frequency, and they
# must come out with OPPOSITE signs -- the source's "A" at 0.80 can only be the
# physical allele the target calls T.  A label-comparing implementation calls
# both EXACT.  d_at_exact_sw (effect T, freq 0.20) is the mirror: the labels
# say "reversed", the frequencies say the source named the target's own A1, and
# the frequencies are right.
#
# The b_* rows pin the MAF threshold on BOTH sides.  A target frequency of
# exactly 0.5 makes the two distances equal and used to fall through to
# RESOLVED_REVERSED -- negating a beta on zero evidence.  0.4 and 0.6 are the
# same MAF and must decide the same way; the frequency parser is documented as
# "ok to lose last ~2 bits of precision", so that is not free.
#
# Verified by mutation: fails against a build with the frequency comparison
# inverted (d_at_exact_lo/d_at_rev_lo swap), and against a build that ignores
# frequencies and compares labels (d_at_rev_lo and d_at_exact_sw invert).
#
# Columns: ID  BIM_A1 BIM_A2  SRC_EFFECT SRC_OTHER  SRC_FREQ  TGT_A1_FREQ
#          WANT_STATUS  WANT_OUT_BETA          (position is 1000 + line number)
# ---------------------------------------------------------------------------
cat > "$WDIR/disc_table.txt" <<'EOF'
d_at_exact_lo	A	T	A	T	0.20	0.20	RESOLVED_EXACT	0.5
d_at_rev_lo	A	T	A	T	0.80	0.20	RESOLVED_REVERSED	-0.5
d_at_exact_sw	A	T	T	A	0.20	0.20	RESOLVED_EXACT	0.5
d_at_rev_sw	A	T	T	A	0.80	0.20	RESOLVED_REVERSED	-0.5
d_cg_exact_lo	C	G	C	G	0.20	0.20	RESOLVED_EXACT	0.5
d_cg_rev_lo	C	G	C	G	0.80	0.20	RESOLVED_REVERSED	-0.5
d_cg_exact_sw	C	G	G	C	0.20	0.20	RESOLVED_EXACT	0.5
d_cg_rev_sw	C	G	G	C	0.80	0.20	RESOLVED_REVERSED	-0.5
d_ctl_nonpal	A	G	A	G	0.80	0.20	EXACT	0.5
b_tgt_tie	A	T	A	T	0.01	0.5	AMBIGUOUS_UNRESOLVED	NA
b_tgt_at_hi	A	T	A	T	0.01	0.6	AMBIGUOUS_UNRESOLVED	NA
b_tgt_at_lo	A	T	A	T	0.01	0.4	AMBIGUOUS_UNRESOLVED	NA
b_tgt_in_hi	A	T	A	T	0.01	0.61	RESOLVED_REVERSED	-0.5
b_tgt_in_lo	A	T	A	T	0.01	0.39	RESOLVED_EXACT	0.5
b_src_at_hi	A	T	A	T	0.6	0.05	AMBIGUOUS_UNRESOLVED	NA
b_src_at_lo	A	T	A	T	0.4	0.05	AMBIGUOUS_UNRESOLVED	NA
b_src_in_hi	A	T	A	T	0.61	0.05	RESOLVED_REVERSED	-0.5
b_src_in_lo	A	T	A	T	0.39	0.05	RESOLVED_EXACT	0.5
EOF

# The .afreq is written with REF = the .bim's A2 and ALT = its A1, so
# ALT_FREQS is the frequency of the .bim's A1 verbatim -- the number the
# resolution actually compares against.
awk -F'\t' -v OFS='\t' '{print "1", $1, "0", 1000+NR, $2, $3}' \
  "$WDIR/disc_table.txt" > "$WDIR/disc.bim"
awk -F'\t' -v OFS='\t' 'BEGIN{print "#CHROM","ID","REF","ALT","ALT_FREQS","OBS_CT"}
                        {print "1", $1, $3, $2, $7, 1000}' \
  "$WDIR/disc_table.txt" > "$WDIR/disc.afreq"
awk -F'\t' -v OFS='\t' '{print "1", 1000+NR, $4, $5, 0.5, $6}' \
  "$WDIR/disc_table.txt" > "$WDIR/disc_src.txt"

bin/adjutant \
  --bim "$WDIR/disc.bim" \
  --pgs-model "$WDIR/disc_src.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/disc.afreq" --resolve-ambiguous \
  --out "$WDIR/disc" > "$WDIR/disc_cli.txt" 2>&1

C_D_POS=$(col "$WDIR/disc.match_report" POS)
C_D_STATUS=$(col "$WDIR/disc.match_report" STATUS)
C_D_BETA=$(col "$WDIR/disc.match_report" OUT_BETA)

# Case 41n: every row's status AND emitted beta, against the table above.
DISC_BAD=$(awk -F'\t' -v p="$C_D_POS" -v st="$C_D_STATUS" -v b="$C_D_BETA" '
  NR==FNR { k=1000+FNR; id[k]=$1; w_st[k]=$8; w_b[k]=$9; n_want++; next }
  FNR>1 {
    ++n_got
    k=$p
    if (!(k in id)) { print "row at POS "k" is not in the table"; next }
    if ($st != w_st[k] || $b != w_b[k])
      print id[k]": got "$st" "$b"  want "w_st[k]" "w_b[k]
  }
  END { if (n_got != n_want) print "report has "n_got" rows, table has "n_want }
  ' "$WDIR/disc_table.txt" "$WDIR/disc.match_report")
if [ -z "$DISC_BAD" ]; then
  echo "PASS  --resolve-ambiguous: 18-row discriminating fixture, every status and beta"
else
  echo "FAIL  --resolve-ambiguous: discriminating fixture:"
  echo "$DISC_BAD" | head -6
  exit 1
fi

# Case 41o: the two rows that carry identical allele labels really do come out
# with opposite signs.  Stated separately from the table so a future edit that
# quietly makes them agree is caught by name rather than by row count.
DISC_LO=$(awk -F'\t' -v p="$C_D_POS" -v b="$C_D_BETA" '$p==1001{print $b}' "$WDIR/disc.match_report")
DISC_HI=$(awk -F'\t' -v p="$C_D_POS" -v b="$C_D_BETA" '$p==1002{print $b}' "$WDIR/disc.match_report")
if [ "$DISC_LO" = "0.5" ] && [ "$DISC_HI" = "-0.5" ]; then
  echo "PASS  --resolve-ambiguous: same alleles, source freq 0.20 vs 0.80 -> +0.5 vs -0.5"
else
  echo "FAIL  --resolve-ambiguous: identical-label rows gave betas '$DISC_LO' and '$DISC_HI', want 0.5 and -0.5"
  exit 1
fi

# Case 41p: a genuinely duplicated bi-allelic .afreq ID.  The FIRST row wins,
# so extra rows for an ID already present cannot change a single answer.
#
# This went untested because Tests/score_fixtures/freqs.afreq's 9 repeated IDs
# are each one bi-allelic plus one multi-allelic row, and multi-allelic rows
# are dropped at parse time -- so no duplicate was ever stored.  Verified by
# mutation: inserting the .afreq IDs in descending order, which is what the
# code used to do, hands the win to the LAST row and turns d_at_exact_lo's
# +0.5 into a drop.
awk -F'\t' -v OFS='\t' '
  { print }
  $2=="d_at_exact_lo" && !done { $5="0.80"; print; $5="0.5"; print; done=1 }
  ' "$WDIR/disc.afreq" > "$WDIR/disc_dup.afreq"
bin/adjutant \
  --bim "$WDIR/disc.bim" \
  --pgs-model "$WDIR/disc_src.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/disc_dup.afreq" --resolve-ambiguous \
  --out "$WDIR/disc_dup" > "$WDIR/disc_dup_cli.txt" 2>&1
if diff -q "$WDIR/disc.match_report" "$WDIR/disc_dup.match_report" > /dev/null; then
  echo "PASS  --resolve-ambiguous: duplicate bi-allelic .afreq IDs keep the first row"
else
  echo "FAIL  --resolve-ambiguous: a duplicate .afreq ID changed the answer:"
  diff "$WDIR/disc.match_report" "$WDIR/disc_dup.match_report" | head -6
  exit 1
fi

# Case 41q: a .bim and an .afreq on different ID schemes join to nothing.  The
# only symptom used to be "0 resolved", which reads like "no palindromes".
awk -F'\t' -v OFS='\t' 'NR==1{print;next}{$2="rs"NR; print}' \
  "$WDIR/disc.afreq" > "$WDIR/disc_noid.afreq"
bin/adjutant \
  --bim "$WDIR/disc.bim" \
  --pgs-model "$WDIR/disc_src.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/disc_noid.afreq" --resolve-ambiguous \
  --out "$WDIR/disc_noid" > "$WDIR/disc_noid_cli.txt" 2>&1
NOID_RES=$(awk -v st="$C_D_STATUS" '$st ~ /^RESOLVED_/' "$WDIR/disc_noid.match_report" | wc -l | tr -d ' ')
if saidflat "$WDIR/disc_noid_cli.txt" "none of the .bim variant IDs appear" && [ "$NOID_RES" -eq 0 ]; then
  echo "PASS  --resolve-ambiguous: an ID join that matches nothing warns"
else
  echo "FAIL  --resolve-ambiguous: silent on an ID join that matched nothing ($NOID_RES resolved)"
  grep -i "warn" "$WDIR/disc_noid_cli.txt"; exit 1
fi

# Case 41r: a .afreq where only SOME rows are malformed used to be accepted in
# silence -- only an all-bad file tripped "no usable rows".
awk -F'\t' -v OFS='\t' 'NR==2{print $1,$2;next} NR==3{$5="not-a-number"} {print}' \
  "$WDIR/disc.afreq" > "$WDIR/disc_partial.afreq"
bin/adjutant \
  --bim "$WDIR/disc.bim" \
  --pgs-model "$WDIR/disc_src.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/disc_partial.afreq" --resolve-ambiguous \
  --out "$WDIR/disc_partial" > "$WDIR/disc_partial_cli.txt" 2>&1
# The filename in the warning is word-wrapped onto its own line by logerrprintfww,
# so match the count, which is what this case is about.
if saidflat "$WDIR/disc_partial_cli.txt" "2 of 18 rows"; then
  echo "PASS  --resolve-ambiguous: malformed .afreq rows are counted, not swallowed"
else
  echo "FAIL  --resolve-ambiguous: no skipped-row count for a partially malformed .afreq"
  grep -i "warn" "$WDIR/disc_partial_cli.txt"; exit 1
fi

# Case 41s: --read-freq is dead weight without --resolve-ambiguous -- the file is
# never opened.  It must say so, and the output must be provably unchanged.
#
# The warning must scope itself to --pgs-model, which is the only mode it can
# fire in.  --score genuinely does read --read-freq, to mean-impute missing
# genotypes, so an unqualified "--read-freq is ignored" invites exactly the
# wrong generalisation from a user who reads it once.
bin/adjutant \
  --bim "$WDIR/disc.bim" \
  --pgs-model "$WDIR/disc_src.txt" 1 2 3 4 5 6 \
  --out "$WDIR/disc_nofreq" > "$WDIR/disc_nofreq_cli.txt" 2>&1
bin/adjutant \
  --bim "$WDIR/disc.bim" \
  --pgs-model "$WDIR/disc_src.txt" 1 2 3 4 5 6 \
  --read-freq "$WDIR/disc_table.txt" \
  --out "$WDIR/disc_junkfreq" > "$WDIR/disc_junkfreq_cli.txt" 2>&1
# Wording follows the code: the guard now covers --pgs-models-list too, so it
# names the mode rather than the single flag.
if saidflat "$WDIR/disc_junkfreq_cli.txt" "in harmonisation mode --read-freq is ignored without --resolve-ambiguous" &&
   diff -q "$WDIR/disc_nofreq.match_report" "$WDIR/disc_junkfreq.match_report" > /dev/null; then
  echo "PASS  --resolve-ambiguous: --read-freq without it warns, and is genuinely unused"
else
  echo "FAIL  --resolve-ambiguous: --read-freq without it was silent or changed the output"
  grep -i "warn" "$WDIR/disc_junkfreq_cli.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 42-54: --folds, k-fold cross-validated R2.
#
# The load-bearing case is 42: k = 1 must reproduce the pre-cross-validation
# output byte for byte, which is what makes every number in every case above
# still mean what it did.  The rest check that k > 1 changes the R2 family and
# only the R2 family, that the interval keeps containing its point estimate,
# and that a k the data cannot support is refused before any scoring happens.
# ---------------------------------------------------------------------------
C_FOLDS=$(col "$BEST" FOLDS)

# Case 42: --folds 1 is byte-identical to no --folds at all, on QT and on CC
# with --prevalence (which exercises the liability columns too).
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --folds 1 --out "$WDIR/qt_f1" > /dev/null 2>&1

bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno_cc.txt" --pheno-name CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --prevalence 0.05 --folds 1 --out "$WDIR/cc_f1" > /dev/null 2>&1

if cmp -s "$BEST" "$WDIR/qt_f1.best_pgs" && \
   cmp -s "$WDIR/gb_cc_prev.best_pgs" "$WDIR/cc_f1.best_pgs"; then
  echo "PASS  --folds: k=1 is byte-identical to no --folds (QT and CC+prevalence)"
else
  echo "FAIL  --folds: k=1 differs from no --folds:"
  diff "$BEST" "$WDIR/qt_f1.best_pgs" | head -6
  diff "$WDIR/gb_cc_prev.best_pgs" "$WDIR/cc_f1.best_pgs" | head -6
  exit 1
fi

# The cross-validated runs the next several cases read.
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --folds 5 --out "$WDIR/qt_f5" > /dev/null 2>&1

bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/pheno_cc.txt" --pheno-name CC_PHENO \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --prevalence 0.05 --folds 5 --out "$WDIR/cc_f5" > /dev/null 2>&1

# Case 43: FOLDS reports the k that was used, on both trait types.
BAD_FOLDS=$(cat "$WDIR/qt_f5.best_pgs" "$WDIR/cc_f5.best_pgs" | grep -v '^#' | \
  awk -v c="$C_FOLDS" '$c != "5" {print NR": "$c}')
if [ -z "$BAD_FOLDS" ]; then
  echo "PASS  --folds: FOLDS column reports the k used"
else
  echo "FAIL  --folds: FOLDS column should read 5:"; printf '%s\n' "$BAD_FOLDS"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 44: containment under cross-validation, and the NA rule that makes it
# possible.
#
# Out of fold R2_INCR can legitimately come out negative -- that is the honest
# answer that the score predicts held-out samples worse than nothing.  The
# interval is on a correlation and its square is non-negative, so there is no
# interval to report there; both bounds must be NA rather than a positive
# interval sitting entirely above its own point estimate.  Where R2_INCR is
# non-negative the interval must contain it, exactly as at k=1.
# ---------------------------------------------------------------------------
CV_BAD_CI=$(cat "$WDIR/qt_f5.best_pgs" "$WDIR/cc_f5.best_pgs" | grep -v '^#' | \
  awk -v id="$C_PGS_ID" -v l="$C_R2_INCR_L95" -v v="$C_R2_INCR" -v u="$C_R2_INCR_U95" '
    $v+0 < 0 {
      if ($l != "NA" || $u != "NA") print $id": negative R2_INCR "$v" with bounds ["$l", "$u"]"
      next
    }
    $l=="NA" || $u=="NA"            { print $id": NA bound on non-negative R2_INCR "$v; next }
    !($l+0 <= $v+0 && $v+0 <= $u+0) { print $id": ["$l", "$v", "$u"]" }')
if [ -z "$CV_BAD_CI" ]; then
  echo "PASS  --folds: CI contains R2_INCR, and negative R2_INCR gets NA bounds"
else
  echo "FAIL  --folds: CI/point-estimate disagreement under cross-validation:"
  printf '%s\n' "$CV_BAD_CI"; exit 1
fi

# Case 44b: the liability pair follows the observed-scale pair, NA and all.
CV_BAD_R2L=$(grep -v '^#' "$WDIR/cc_f5.best_pgs" | \
  awk -v id="$C_PGS_ID" -v ov="$C_R2_INCR" \
      -v l="$C_R2L_INCR_L95" -v v="$C_R2L_INCR" -v u="$C_R2L_INCR_U95" '
    $ov+0 < 0 { if ($l != "NA" || $u != "NA") print $id": R2L bounds not NA"; next }
    !($l+0 <= $v+0 && $v+0 <= $u+0) { print $id": ["$l", "$v", "$u"]" }')
if [ -z "$CV_BAD_R2L" ]; then
  echo "PASS  --folds: R2L_INCR interval tracks the observed-scale one"
else
  echo "FAIL  --folds: R2L_INCR interval is wrong under cross-validation:"
  printf '%s\n' "$CV_BAD_R2L"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 45: R2_INCR is still R2_FULL - R2_NULL under cross-validation.
#
# Case 21b asserts this at k=1.  It has to survive k>1 as well, and it only can
# if all three come from the same pooled out-of-fold set against the same TSS.
# If this ever fails, the three were computed over different sample sets --
# which is a bug, not a reason to relax the tolerance.
# ---------------------------------------------------------------------------
CV_BAD_INCR=$(cat "$WDIR/qt_f5.best_pgs" "$WDIR/cc_f5.best_pgs" | grep -v '^#' | \
  awk -v id="$C_PGS_ID" -v t="$C_TYPE" \
      -v rn="$C_R2_NULL" -v rf="$C_R2_FULL" -v ri="$C_R2_INCR" \
      -v an="$C_AUC_NULL" -v af="$C_AUC_FULL" -v ai="$C_AUC_INCR" '
    function abs(x) { return x < 0 ? -x : x }
    function off(lo, hi, d,    e) {
      e = abs((hi - lo) - d)
      return e > 1e-4 * (abs(lo) + abs(hi) + abs(d))
    }
    off($rn, $rf, $ri) { print $id" R2: "$rf" - "$rn" != "$ri }
    $t=="CC" && off($an, $af, $ai) { print $id" AUC: "$af" - "$an" != "$ai }')
if [ -z "$CV_BAD_INCR" ]; then
  echo "PASS  --folds: R2_INCR and AUC_INCR are still the differences they claim"
else
  echo "FAIL  --folds: an _INCR column is not the difference of its pair:"
  printf '%s\n' "$CV_BAD_INCR"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 46: on the causal model, cross-validation moves the R2 family and leaves
# the coefficient alone.  Both fits explain less of the data they did not see:
# that is what is asserted, per model.  R2_INCR itself is a difference of two
# quantities that both shrink, so its direction is not fixed for a single
# model -- the shrinkage that matters is the winner-of-N one, and that is
# Case 48.
# ---------------------------------------------------------------------------
read W2_IN_FULL W2_IN_NULL W2_IN_BETA W2_IN_AUC <<EOF
$(awk -F'\t' -v id="$C_PGS_ID" -v f="$C_R2_FULL" -v n="$C_R2_NULL" -v b="$C_BETA" -v a="$C_AUC_FULL" \
   '$id=="W2"{print $f, $n, $b, $a}' "$WDIR/gb_cc_prev.best_pgs")
EOF
read W2_CV_FULL W2_CV_NULL W2_CV_BETA W2_CV_AUC W2_CV_INCR <<EOF
$(awk -F'\t' -v id="$C_PGS_ID" -v f="$C_R2_FULL" -v n="$C_R2_NULL" -v b="$C_BETA" -v a="$C_AUC_FULL" -v i="$C_R2_INCR" \
   '$id=="W2"{print $f, $n, $b, $a, $i}' "$WDIR/cc_f5.best_pgs")
EOF
CV_STRONG=$(awk -v inf="$W2_IN_FULL" -v cvf="$W2_CV_FULL" \
                -v inn="$W2_IN_NULL" -v cvn="$W2_CV_NULL" -v cvi="$W2_CV_INCR" \
  'BEGIN{
     if (!(cvf+0 < inf+0)) print "R2_FULL did not shrink: "cvf" >= "inf
     if (!(cvn+0 < inn+0)) print "R2_NULL did not shrink: "cvn" >= "inn
     if (!(cvi+0 > 0.2))   print "R2_INCR no longer clearly positive: "cvi
   }')
if [ -z "$CV_STRONG" ]; then
  echo "PASS  --folds: on the causal model both R2s shrink and R2_INCR stays large"
else
  echo "FAIL  --folds: cross-validation on the causal model:"
  printf '%s\n' "$CV_STRONG"; exit 1
fi

# Case 47: BETA and the AUC family are documented as in-sample, so --folds must
# not move them.  AUC comes from the logistic fit, which is not refit per fold;
# reporting a cross-validated R2 beside an in-sample AUC is only defensible
# because it is stated, so it is also asserted.
if [ "$W2_IN_BETA" = "$W2_CV_BETA" ] && [ "$W2_IN_AUC" = "$W2_CV_AUC" ]; then
  echo "PASS  --folds: BETA and AUC_FULL are unchanged by --folds (in-sample by design)"
else
  echo "FAIL  --folds: BETA/AUC moved: beta $W2_IN_BETA -> $W2_CV_BETA, auc $W2_IN_AUC -> $W2_CV_AUC"
  exit 1
fi

# ---------------------------------------------------------------------------
# Case 48: 50 models of pure noise against an unrelated phenotype.  This is the
# case the whole flag exists for: ranking N models and reporting the winner's
# in-sample R2 overstates it, and the overstatement grows with N.  The winner
# still wins -- cross-validation fixes the estimate, not the selection -- but
# its R2_INCR must come down, and most of the field must fall below zero, which
# is the answer an in-sample R2 cannot give.
# ---------------------------------------------------------------------------
awk 'BEGIN{srand(1234); h="ID\tA1"; for(k=1;k<=50;k++) h=h"\tN"k; print h}
     !/^#/{line=$3"\t"$5; for(k=1;k<=50;k++) line=line"\t"(rand()-0.5); print line}' \
  "$PFILE_PFX.pvar" | head -400 > "$WDIR/noise50.txt"

for kk in 1 5; do
  bin/adjutant \
    --pfile "$PFILE_PFX" vzs \
    --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
    --covar "$WDIR/covar.txt" \
    --score "$WDIR/noise50.txt" 1 2 header-read --score-col-nums 3-52 \
    --give-best-pgs --folds $kk --out "$WDIR/noise50_k$kk" > /dev/null 2>&1
done

NOISE_WINNER=$(awk -F'\t' -v c="$C_PGS_ID" 'NR==2{print $c}' "$WDIR/noise50_k1.best_pgs")
NOISE_IN=$(awk -F'\t' -v id="$C_PGS_ID" -v v="$C_R2_INCR" -v w="$NOISE_WINNER" \
           '$id==w{print $v}' "$WDIR/noise50_k1.best_pgs")
NOISE_CV=$(awk -F'\t' -v id="$C_PGS_ID" -v v="$C_R2_INCR" -v w="$NOISE_WINNER" \
           '$id==w{print $v}' "$WDIR/noise50_k5.best_pgs")
NOISE_NEG=$(awk -F'\t' -v v="$C_R2_INCR" 'NR>1 && $v+0 < 0' "$WDIR/noise50_k5.best_pgs" | wc -l | tr -d ' ')
NOISE_BAD=$(awk -v i="$NOISE_IN" -v c="$NOISE_CV" -v n="$NOISE_NEG" \
  'BEGIN{
     if (!(c+0 < i+0)) print "winner R2_INCR did not shrink: "c" >= "i
     if (!(n+0 >= 10)) print "only "n" of 50 noise models fall below zero out of fold"
   }')
if [ -z "$NOISE_BAD" ]; then
  echo "PASS  --folds: 50-noise-model winner $NOISE_WINNER shrinks $NOISE_IN -> $NOISE_CV ($NOISE_NEG/50 negative)"
else
  echo "FAIL  --folds: 50-noise-model case:"; printf '%s\n' "$NOISE_BAD"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 49: an imbalanced binary phenotype (39 cases in 2504) with --folds 10.
# Folds are stratified on case/control, so every fold still holds about four
# cases and the null model can be fitted on all ten training sets.  Without the
# stratification an unlucky split leaves a fold with none, and the run either
# fails or returns a meaningless R2.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tCC_IMB"} NR>1{print $1"\t"((NR%64==0)?2:1)}' \
  "$PSAM" > "$WDIR/pheno_imb.txt"

if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/pheno_imb.txt" --pheno-name CC_IMB \
     --covar "$WDIR/covar.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --folds 10 --out "$WDIR/imb_f10" > "$WDIR/imb_f10_out.txt" 2>&1
then
  IMB_BAD=$(grep -v '^#' "$WDIR/imb_f10.best_pgs" | \
    awk -v id="$C_PGS_ID" -v b="$C_BETA" -v nc="$C_N_CASE" -v f="$C_FOLDS" \
        '$b=="NA" || $nc!="39" || $f!="10" {print $id": BETA="$b" N_CASE="$nc" FOLDS="$f}')
  IMB_ROWS=$(grep -v '^#' "$WDIR/imb_f10.best_pgs" | wc -l | tr -d ' ')
  if [ -z "$IMB_BAD" ] && [ "$IMB_ROWS" -eq 3 ]; then
    echo "PASS  --folds: imbalanced CC (39 cases) fits in all 10 folds"
  else
    echo "FAIL  --folds: imbalanced CC run produced $IMB_ROWS rows:"
    printf '%s\n' "$IMB_BAD"; exit 1
  fi
else
  echo "FAIL  --folds: imbalanced CC with --folds 10 aborted"
  cat "$WDIR/imb_f10_out.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 50: --folds 0 is refused at parse time -- before the pgen, phenotypes and
# covariates are even loaded, never mind scored.  k = 0 is wrong on its face,
# unlike a k that is only too large for this dataset, so it belongs with the
# other static argument checks; "no samples loaded" is what distinguishes the
# two, and the absent .sscore keeps the pre-flight promise either way.
# ---------------------------------------------------------------------------
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --covar "$WDIR/covar.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --folds 0 --out "$WDIR/f0" > "$WDIR/f0_out.txt" 2>&1
then
  echo "FAIL  --folds 0: should have been refused"; exit 1
else
  if grep -qi "folds" "$WDIR/f0_out.txt" && [ ! -e "$WDIR/f0.sscore" ] &&
     ! saidflat "$WDIR/f0_out.txt" "samples ("; then
    echo "PASS  --folds 0: refused at parse time, before anything is loaded"
  else
    echo "FAIL  --folds 0: refused too late or with the wrong message"
    cat "$WDIR/f0_out.txt"; exit 1
  fi
fi

# Case 51: --folds greater than half the complete-case count is refused, also
# before scoring.  2504 complete cases, so 1253 is one past the limit.
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --covar "$WDIR/covar.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --folds 1253 --out "$WDIR/fbig" > "$WDIR/fbig_out.txt" 2>&1
then
  echo "FAIL  --folds 1253: should have been refused (2504 complete cases)"; exit 1
else
  if grep -qi "folds" "$WDIR/fbig_out.txt" && [ ! -e "$WDIR/fbig.sscore" ]; then
    echo "PASS  --folds: k above half the sample count is refused before scoring"
  else
    echo "FAIL  --folds 1253: refused too late or with the wrong message"
    cat "$WDIR/fbig_out.txt"; exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Case 51b: a k that clears the fold-size floor but starves the training set is
# refused too, and also before scoring.  Eight complete cases over two folds
# leaves four training samples for three predictors, which is one short of what
# PgsEvalOlsSetFolds needs -- while the floor Case 51 tests waves it through
# (2 <= 8 / 2).  Before the pre-flight knew the real constraint this run wrote a
# whole .sscore and only then failed inside the fitter, which is precisely the
# contract BestPgsValidate exists to keep, so the absent .sscore is the half
# that matters.  'no-mean-imputation' is only there to get past the 50-sample
# allele-frequency floor, which would otherwise stop the run for an unrelated
# reason and make the .sscore assertion vacuous.
# ---------------------------------------------------------------------------
awk 'NR>1 && NR<=9 {print $1}' "$PSAM" > "$WDIR/keep8.txt"
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs --keep "$WDIR/keep8.txt" \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --covar "$WDIR/covar.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read no-mean-imputation --score-col-nums 3-5 \
     --give-best-pgs --folds 2 --out "$WDIR/fsmall" > "$WDIR/fsmall_out.txt" 2>&1
then
  echo "FAIL  --folds 2 on 8 complete cases: should have been refused"; exit 1
else
  # Unwrapped, because logerrprintfww breaks the message across lines.
  #
  # The tail is pinned as well as the head.  The message used to end "too few to
  # fit 4 predictors on ... or drop a covariate", and both halves were wrong: a
  # constant covariate is dropped by the tool itself, so the fit that would run
  # here uses fewer predictors than the check counts, and advising the user to
  # drop one is advice about work already done.  The counts below are the two
  # the message may now state -- 2 covariates plus the intercept and the PGS,
  # needing 5 training samples -- and both are arithmetic on this fixture, so a
  # message that drifts back to naming the predictor count fails here.
  FSMALL_MSG=$(tr '\n' ' ' < "$WDIR/fsmall_out.txt")
  FSMALL_OK=0
  case "$FSMALL_MSG" in
    *"--folds 2 leaves only 4 samples in the smallest training set"*"intercept, 2 covariate(s) and the PGS there needs at least 5"*)
      case "$FSMALL_MSG" in
        *"drop a covariate"*) FSMALL_OK=0 ;;
        *)                    FSMALL_OK=1 ;;
      esac ;;
  esac
  if [ "$FSMALL_OK" = 1 ] && [ ! -e "$WDIR/fsmall.sscore" ]; then
    echo "PASS  --folds: a k that starves the training set is refused before scoring"
  else
    echo "FAIL  --folds 2 on 8 complete cases: refused too late or with the wrong message"
    cat "$WDIR/fsmall_out.txt"; exit 1
  fi
fi

# Case 52: --folds without --give-best-pgs is refused; there is nothing for it
# to cross-validate.
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --folds 5 --out "$WDIR/fng" > "$WDIR/fng_out.txt" 2>&1
then
  echo "FAIL  --folds without --give-best-pgs: should have been refused"; exit 1
else
  if grep -qi "give-best-pgs" "$WDIR/fng_out.txt"; then
    echo "PASS  --folds: refused without --give-best-pgs"
  else
    echo "FAIL  --folds: refused but wrong error message"
    cat "$WDIR/fng_out.txt"; exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Case 53: folds are seeded, so a rerun of the same command reproduces the same
# numbers and a different --seed produces different ones.  Both halves are
# needed: identical output alone would also be what a fold assignment that
# ignored the seed entirely produces.
# ---------------------------------------------------------------------------
for tag in a b c; do
  case $tag in
    a) SEEDARG="" ;;
    b) SEEDARG="" ;;
    c) SEEDARG="--seed 987654321" ;;
  esac
  bin/adjutant \
    --pfile "$PFILE_PFX" vzs \
    --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
    --covar "$WDIR/covar.txt" \
    --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
    --give-best-pgs --folds 5 $SEEDARG --out "$WDIR/seed_$tag" > /dev/null 2>&1
done
if cmp -s "$WDIR/seed_a.best_pgs" "$WDIR/seed_b.best_pgs" && \
   ! cmp -s "$WDIR/seed_a.best_pgs" "$WDIR/seed_c.best_pgs"; then
  echo "PASS  --folds: fold assignment is reproducible and --seed changes it"
else
  echo "FAIL  --folds: reproducibility or seed sensitivity is broken"
  diff "$WDIR/seed_a.best_pgs" "$WDIR/seed_b.best_pgs" | head -4
  exit 1
fi

# Case 54: --folds works on the --sscore input path too, which reaches the same
# writer by a different route and does its own validation.
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --folds 5 --out "$WDIR/viafile_f5" > /dev/null 2>&1

grep -v '^#' "$WDIR/qt_f5.best_pgs" > "$WDIR/qt_f5_data.txt"
grep -v '^#' "$WDIR/viafile_f5.best_pgs" > "$WDIR/viafile_f5_data.txt"
if diff "$WDIR/qt_f5_data.txt" "$WDIR/viafile_f5_data.txt" > /dev/null 2>&1; then
  echo "PASS  --folds: --sscore and --pfile routes agree under cross-validation"
else
  echo "FAIL  --folds: --sscore and --pfile routes differ under cross-validation:"
  diff "$WDIR/qt_f5_data.txt" "$WDIR/viafile_f5_data.txt" | head -6
  exit 1
fi

# Case 54b: --folds above half the joined sample count is refused on the
# --sscore path as well, which does not go through BestPgsValidate.
if bin/adjutant \
     --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
     --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
     --covar "$WDIR/covar.txt" \
     --give-best-pgs --folds 1253 --out "$WDIR/viafile_fbig" > "$WDIR/viafile_fbig_out.txt" 2>&1
then
  echo "FAIL  --sscore --folds 1253: should have been refused"; exit 1
else
  if grep -qi "folds" "$WDIR/viafile_fbig_out.txt"; then
    echo "PASS  --folds: --sscore path refuses a k the joined sample set cannot support"
  else
    echo "FAIL  --sscore --folds: refused but wrong error message"
    cat "$WDIR/viafile_fbig_out.txt"; exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Case 55: --folds does not turn a would-be-skipped phenotype into a fatal
# error.
#
# BADCC has 6 non-missing values and all of them are controls, so it is skipped
# with a warning and the run continues on QT1.  The fold check used to sit above
# that skip, so --folds 5 made the same phenotype abort the whole run with
# "needs at least 10 complete-case observations" -- a k that the phenotype was
# never going to reach the fitter to use.  Both halves are asserted: the run
# survives, and QT1 alone is in the output.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tQT1\tBADCC"}
     NR>1{print $1"\t"(NR%97)/13.0"\t"((NR<=7)?1:"NA")}' "$PSAM" > "$WDIR/skipcc.txt"
if bin/adjutant \
     --pfile "$PFILE_PFX" vzs \
     --pheno "$WDIR/skipcc.txt" --pheno-name QT1 BADCC \
     --covar "$WDIR/covar.txt" \
     --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
     --give-best-pgs --folds 5 --out "$WDIR/skipcc" > "$WDIR/skipcc_out.txt" 2>&1
then
  SKIPCC_PHENOS=$(grep -v '^#' "$WDIR/skipcc.best_pgs" | awk -v c="$C_PHENO" '{print $c}' | sort -u | tr '\n' ' ')
  if [ "$SKIPCC_PHENOS" = "QT1 " ] && saidflat "$WDIR/skipcc_out.txt" "no cases in the complete-case"; then
    echo "PASS  --folds: an unusable phenotype is still skipped, not fatal"
  else
    echo "FAIL  --folds: expected QT1 rows only and a 'no cases' warning, got '$SKIPCC_PHENOS'"
    cat "$WDIR/skipcc_out.txt"; exit 1
  fi
else
  echo "FAIL  --folds: a phenotype that would be skipped aborted the run"
  cat "$WDIR/skipcc_out.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 56: --sscore refuses a constant phenotype, as --pfile already did.
#
# Without this the fitters accept it: PgsEvalOlsInit's `tss > 0` guards make
# R2_NULL, R2_FULL and R2_INCR all exactly 0, and the CI columns then print a
# Fisher z half-width at r = 0 -- a bound on an R2 that does not exist.  A
# fabricated interval is worse than no row, so the row must be gone entirely.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tQT1\tCONSTQT"} NR>1{print $1"\t"(NR%97)/13.0"\t7"}' \
  "$PSAM" > "$WDIR/constqt.txt"
if bin/adjutant \
     --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
     --pheno "$WDIR/constqt.txt" --pheno-name QT1 CONSTQT \
     --covar "$WDIR/covar.txt" \
     --give-best-pgs --out "$WDIR/ss_const" > "$WDIR/ss_const_out.txt" 2>&1
then
  CONST_PHENOS=$(grep -v '^#' "$WDIR/ss_const.best_pgs" | awk -v c="$C_PHENO" '{print $c}' | sort -u | tr '\n' ' ')
  if [ "$CONST_PHENOS" = "QT1 " ] && saidflat "$WDIR/ss_const_out.txt" "'CONSTQT' is constant; skipping"; then
    echo "PASS  --sscore: a constant phenotype is skipped, not given a CI"
  else
    echo "FAIL  --sscore: expected QT1 rows only and a constant-phenotype warning, got '$CONST_PHENOS'"
    cat "$WDIR/ss_const_out.txt"; exit 1
  fi
else
  echo "FAIL  --sscore: the constant-phenotype run aborted"
  cat "$WDIR/ss_const_out.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 57-58: --sscore treats a missing phenotype value as missing.
#
# It used to coerce one to a control -- `val = (val > 1.5) ? 1.0 : 0.0` turns
# NaN into 0 -- so N counted samples that had no phenotype, N_CASE described a
# sample set that did not exist, and a missing quantitative value poisoned its
# whole column with NaN.  --folds made it load-bearing: both the pre-flight k
# check and the stratified fold assignment ran on the fabricated labels.
#
# The check that means the most is agreement with the --pfile route, which has
# always done the complete-case selection: same phenotype file, same covariates,
# same scores, so the two routes must now produce the same rows to the byte.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tMISSCC\tMISSQT"}
     NR>1{miss=(NR%11==0); print $1"\t"(miss?"NA":((NR%3==0)?2:1))"\t"(miss?"NA":(NR%97)/13.0)}' \
  "$PSAM" > "$WDIR/missing.txt"

bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --pheno "$WDIR/missing.txt" --pheno-name MISSCC MISSQT \
  --covar "$WDIR/covar.txt" \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --folds 5 --out "$WDIR/miss_direct" > "$WDIR/miss_direct_out.txt" 2>&1

bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/missing.txt" --pheno-name MISSCC MISSQT \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --folds 5 --out "$WDIR/miss_viafile" > "$WDIR/miss_viafile_out.txt" 2>&1

grep -v '^#' "$WDIR/miss_direct.best_pgs"  > "$WDIR/miss_direct_data.txt"
grep -v '^#' "$WDIR/miss_viafile.best_pgs" > "$WDIR/miss_viafile_data.txt"
if diff "$WDIR/miss_direct_data.txt" "$WDIR/miss_viafile_data.txt" > /dev/null 2>&1; then
  echo "PASS  --sscore: missing phenotype values give the same rows as --pfile"
else
  echo "FAIL  --sscore: missing phenotype values still differ from --pfile:"
  diff "$WDIR/miss_direct_data.txt" "$WDIR/miss_viafile_data.txt" | head -6
  exit 1
fi

# Case 58: N is the complete-case count on both routes, and the --sscore route
# says how many it dropped.  2504 joined samples, every 11th missing, so 227
# are excluded and 2277 remain -- asserted as literals because the arithmetic is
# the point.  Without this, "the two routes agree" would also pass if both
# reverted to 2504.
MISS_N=$(grep -v '^#' "$WDIR/miss_viafile.best_pgs" | awk -v c="$C_N" '{print $c}' | sort -u | tr '\n' ' ')
if [ "$MISS_N" = "2277 " ] &&
   saidflat "$WDIR/miss_viafile_out.txt" "227 of the 2504 joined samples"; then
  echo "PASS  --sscore: N is the complete-case count and the drop is reported"
else
  echo "FAIL  --sscore: expected N=2277 on every row and a 227-dropped warning, got '$MISS_N'"
  cat "$WDIR/miss_viafile_out.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 59: a ranking in which nothing ranks says so.
#
# The output's shape -- one group per phenotype, best first -- reads as "this
# one", and a cross-validated R2_INCR can be negative for every model, which
# means none of them predicts out of sample.  A noise phenotype is exactly that
# case: at k = 1 all three in-sample R2_INCRs are positive and no warning is
# due, at k = 5 all three go negative and the top row is merely the least bad.
# The two halves together are what make the warning conditional rather than
# unconditional.
# ---------------------------------------------------------------------------
awk 'BEGIN{srand(3); OFS="\t"; print "#IID","NOISECC"} NR>1{print $1, (rand()<0.35)?2:1}' \
  "$PSAM" > "$WDIR/noise_cc.txt"
for k in 1 5; do
  bin/adjutant \
    --pfile "$PFILE_PFX" vzs \
    --pheno "$WDIR/noise_cc.txt" --pheno-name NOISECC \
    --covar "$WDIR/covar.txt" \
    --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
    --give-best-pgs --folds $k --out "$WDIR/noise_$k" > "$WDIR/noise_${k}_out.txt" 2>&1
done
NEG_CT=$(grep -v '^#' "$WDIR/noise_5.best_pgs" | awk -v c="$C_R2_INCR" '$c < 0' | wc -l | tr -d ' ')
if [ "$NEG_CT" -eq 3 ] &&
   saidflat "$WDIR/noise_5_out.txt" "no model has a positive R2_INCR for phenotype" &&
   ! saidflat "$WDIR/noise_1_out.txt" "no model has a positive R2_INCR"; then
  echo "PASS  --folds: an all-negative ranking is called out, a positive one is not"
else
  echo "FAIL  --folds: expected 3 negative R2_INCRs at k=5 with a warning, got $NEG_CT"
  cat "$WDIR/noise_5_out.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 42-44: --sscore accepts PLINK's "#FID IID" phenotype/covariate header.
#
# PLINK writes phenotype and covariate files two ways -- "#IID <data>" and
# "#FID IID <data>" -- and both are ordinary output of its own tools.  The
# --sscore readers used to assume the first column was always the IID, which
# broke the second form in two ways at once:
#
#   * every column after the first became DATA, so IID itself became a
#     covariate.  Its values are sample names, which scan as missing, so every
#     sample had a "missing covariate" and the run died with "no usable
#     phenotype found" -- a message pointing nowhere near the header.
#   * the join key was the FID, which matches nothing as soon as FID != IID.
#
# The first failure was reported from real use.  The second is quieter: it
# hides completely whenever FID == IID, which is the common case.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#FID\tIID\tQT1"; next} {print $1"\t"$1"\t"$2}' OFS='\t' \
  "$WDIR/pheno.txt" > "$WDIR/pheno_fid.txt"
awk 'NR==1{print "#FID\tIID\tC1\tC2"; next} {print $1"\t"$1"\t"$2"\t"$3}' OFS='\t' \
  "$WDIR/covar.txt" > "$WDIR/covar_fid.txt"
# FID deliberately unlike IID, so a reader that joins on column 1 finds nothing.
awk 'NR==1{print "#FID\tIID\tQT1"; next} {print "FAM"NR"\t"$1"\t"$2}' OFS='\t' \
  "$WDIR/pheno.txt" > "$WDIR/pheno_fid_diff.txt"

bin/adjutant --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pheno_fid.txt" --pheno-name QT1 --covar "$WDIR/covar_fid.txt" \
  --give-best-pgs --out "$WDIR/fidhdr" > "$WDIR/fidhdr_out.txt" 2>&1 || true

# Case 42: the reported failure -- a #FID IID covariate file must not exclude
# every sample.
if [ -s "$WDIR/fidhdr.best_pgs" ] &&
   ! saidflat "$WDIR/fidhdr_out.txt" "no usable phenotype found"; then
  echo "PASS  --sscore: a '#FID IID' covariate file does not exclude every sample"
else
  echo "FAIL  --sscore: '#FID IID' covariate file still excludes everything"
  cat "$WDIR/fidhdr_out.txt"; exit 1
fi

# Case 43: the sharp one.  The header form is presentation, not data, so the
# two spellings of the same table must produce byte-identical results.  A
# reader that mis-assigns columns can still emit rows; only equality catches it.
if diff "$WDIR/fidhdr.best_pgs" "$WDIR/viafile.best_pgs" > /dev/null 2>&1; then
  echo "PASS  --sscore: '#FID IID' and '#IID' headers give identical results"
else
  echo "FAIL  --sscore: header form changed the results:"
  diff "$WDIR/viafile.best_pgs" "$WDIR/fidhdr.best_pgs" | head -6
  exit 1
fi

# Case 44: the join must key on IID, not on the first column.
bin/adjutant --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pheno_fid_diff.txt" --pheno-name QT1 \
  --give-best-pgs --out "$WDIR/fiddiff" > "$WDIR/fiddiff_out.txt" 2>&1 || true
if [ -s "$WDIR/fiddiff.best_pgs" ] &&
   ! saidflat "$WDIR/fiddiff_out.txt" "no samples remain after joining"; then
  echo "PASS  --sscore: joins on IID even when FID differs from it"
else
  echo "FAIL  --sscore: joined on FID rather than IID"
  cat "$WDIR/fiddiff_out.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# The .sscore's OWN "#FID IID" header.
#
# The cases above cover the phenotype and covariate readers.  The .sscore
# reader had the same bug and nothing here caught it, because every .sscore in
# this suite comes from --pfile, whose .psam has no FID -- so the header is
# "#IID ..." and the FID branch is never taken.
#
# --bfile output ALWAYS takes it: a .fam carries FID, so FidColIsRequired is
# true under 'maybefid' and ScoreReport writes the column.  The reader took
# token 0 as the IID regardless, so every row's ID came out as the FID.  With
# --make-bed's default FID of 0 that is 300 samples all called "0", and the run
# died with `duplicate IID '0'`.  adjutant could not read back its own --bfile
# .sscore.
# ---------------------------------------------------------------------------
# Same scores, FID prepended.  All FIDs identical, reproducing --make-bed.
awk -F'\t' 'NR==1{print "#FID\t"substr($0,2); next} {print "0\t"$0}' \
  "$WDIR/direct.sscore" > "$WDIR/ss_fid.sscore"
# And a version where FID differs from IID, so joining on column 1 finds nothing.
awk -F'\t' 'NR==1{print "#FID\t"substr($0,2); next} {print "FAM"NR"\t"$0}' \
  "$WDIR/direct.sscore" > "$WDIR/ss_fiddiff.sscore"

bin/adjutant --sscore "$WDIR/ss_fid.sscore" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/ss_fid" > "$WDIR/ss_fid_out.txt" 2>&1 || true
need "$WDIR/ss_fid.best_pgs" "--sscore: a '#FID IID' .sscore produced nothing"

# The header form is presentation, not data: the rows must be identical to the
# same scores read without an FID column.
grep -v '^#' "$WDIR/viafile.best_pgs" > "$WDIR/ss_nofid_data.txt"
grep -v '^#' "$WDIR/ss_fid.best_pgs"  > "$WDIR/ss_fid_data.txt"
if diff "$WDIR/ss_nofid_data.txt" "$WDIR/ss_fid_data.txt" > /dev/null 2>&1; then
  echo "PASS  --sscore: a '#FID IID' .sscore gives the same rows as '#IID'"
else
  echo "FAIL  --sscore: the .sscore FID column changed the result:"
  diff "$WDIR/ss_nofid_data.txt" "$WDIR/ss_fid_data.txt" | head -6; exit 1
fi

# The sharp one: identical FIDs must not read as duplicate IIDs, and a FID
# unlike the IID must still join.  Both are the failure the bug produced.
bin/adjutant --sscore "$WDIR/ss_fiddiff.sscore" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/ss_fiddiff" > "$WDIR/ss_fiddiff_out.txt" 2>&1 || true
if [ -s "$WDIR/ss_fiddiff.best_pgs" ] &&
   ! saidflat "$WDIR/ss_fid_out.txt" "duplicate IID" &&
   ! saidflat "$WDIR/ss_fiddiff_out.txt" "no samples remain after joining"; then
  echo "PASS  --sscore: .sscore IIDs are read past the FID column, not as it"
else
  echo "FAIL  --sscore: .sscore FID column read as the IID"
  cat "$WDIR/ss_fid_out.txt" "$WDIR/ss_fiddiff_out.txt" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 45-46: --sscore auto-detection must not rank non-score columns.
#
# --score's DEFAULT column set includes 'phenos', so the .sscore it writes is
#
#   #IID  PHENO1  PHENO2  ALLELE_CT  NAMED_ALLELE_DOSAGE_SUM  W1_AVG ...
#
# The old auto-detection scanned left to right for the first name it did not
# recognise and treated EVERYTHING after it as a score.  That assumes the
# non-score columns form a contiguous prefix; 'phenos' sits in the middle, so
# PHENO1 tripped the switch and PHENO1, PHENO2, ALLELE_CT and
# NAMED_ALLELE_DOSAGE_SUM were all ranked as candidate models.
#
# That is not cosmetic.  When the .sscore carries the phenotype under analysis
# -- the normal case when the fileset already holds it -- that column predicts
# itself and wins: measured R2_INCR 0.98 at P = 0, above every real model.
#
# The suite could not have caught this before: 11 of its 12 --sscore cases pass
# --score-col-nums, so the default path was the least covered one.
# ---------------------------------------------------------------------------
bin/adjutant --sscore "$WDIR/direct.sscore" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/autodet" > "$WDIR/autodet_out.txt" 2>&1 || true

# Case 45: no phenotype or bookkeeping column may appear as a model.
AUTO_IDS=$(grep -v '^#' "$WDIR/autodet.best_pgs" 2>/dev/null | cut -f1 | sort -u | tr '\n' ' ')
# `grep -c` exits 1 on zero matches, which is the PASSING outcome here, so
# without the guard set -e kills the suite exactly when the code is correct.
BAD=$(printf '%s\n' $AUTO_IDS | grep -cE '^(PHENO[0-9]*|ALLELE_CT|NAMED_ALLELE_DOSAGE.*|DENOM|SUPERPOP|POPULATION)$' || true)
if [ "$BAD" -eq 0 ]; then
  echo "PASS  --sscore: auto-detection ignores phenotype and bookkeeping columns"
else
  echo "FAIL  --sscore: auto-detection ranked non-score columns: $AUTO_IDS"
  exit 1
fi

# Case 46: the sharp one.  Auto-detection must select exactly the columns the
# explicit range selects -- same models, same numbers.  Counting rows alone
# would pass on a set that is merely the right SIZE.
EXPLICIT_IDS=$(grep -v '^#' "$WDIR/viafile.best_pgs" | cut -f1 | sort -u | tr '\n' ' ')
if [ "$AUTO_IDS" = "$EXPLICIT_IDS" ]; then
  echo "PASS  --sscore: auto-detected columns match --score-col-nums exactly"
else
  echo "FAIL  --sscore: auto='$AUTO_IDS' explicit='$EXPLICIT_IDS'"
  exit 1
fi
# Cases 47-57: --pgs-models-list.
#
# Merging several models into one --score matrix is only meaningful because
# every model is harmonised against the SAME .bim first: two models may name
# opposite effect alleles for one variant, and merging raw files on variant ID
# would invert one column silently.  These cases pin that, the union/zero-fill
# shape, and the guards.
#
# The three models are overlapping SLICES OF ONE SOURCE FILE, which buys the
# sharpest assertion available: where two models carry the same variant they
# carry the same source row, so their merged columns must be equal.  A merge
# that mismatched alleles or mis-indexed the matrix breaks that equality even
# when the row count still looks right.
# ---------------------------------------------------------------------------
MLSRC="$WDIR/pgs_src.txt"
MLN=$(grep -vc '^#' "$MLSRC")
ML1=$((MLN / 3))
ML2=$((MLN * 2 / 3))
awk -v lo=0            -v hi=$ML1 '!/^#/{n++; if(n>lo && n<=hi) print}' "$MLSRC" > "$WDIR/mL_A.txt"
awk -v lo=$((ML1 / 2)) -v hi=$ML2 '!/^#/{n++; if(n>lo && n<=hi) print}' "$MLSRC" > "$WDIR/mL_B.txt"
awk -v lo=$((ML2 / 2)) -v hi=$MLN '!/^#/{n++; if(n>lo && n<=hi) print}' "$MLSRC" > "$WDIR/mL_C.txt"
printf '%s\n# a comment; the indices apply to every listed file\n%s\n%s\n' \
  "$WDIR/mL_A.txt" "$WDIR/mL_B.txt" "$WDIR/mL_C.txt" > "$WDIR/mL_list.txt"
printf '# only a comment, so no paths at all\n' > "$WDIR/mL_empty.txt"

bin/adjutant --bim "$BIM" --pgs-models-list "$WDIR/mL_list.txt" 1 2 3 4 5 \
  --out "$WDIR/mL" > "$WDIR/mL_out.txt" 2>&1

# Case 42: header is ID, A1, then one column per model named from the basename.
MLHDR=$(head -1 "$WDIR/mL.matched")
if [ "$MLHDR" = "$(printf 'ID\tA1\tmL_A\tmL_B\tmL_C')" ]; then
  echo "PASS  --pgs-models-list: header names columns from the file basenames"
else
  echo "FAIL  --pgs-models-list: header is '$MLHDR'"; exit 1
fi

# Case 43: no row is all-zero.  One would mean a variant entered the union
# without any model supplying it.
MLZ=$(awk -F'\t' 'NR>1 && $3==0 && $4==0 && $5==0' "$WDIR/mL.matched" | wc -l | tr -d ' ')
if [ "$MLZ" -eq 0 ]; then
  echo "PASS  --pgs-models-list: no all-zero row in the union"
else
  echo "FAIL  --pgs-models-list: $MLZ all-zero rows"; exit 1
fi

# Case 44: the models genuinely overlap AND genuinely differ, so case 43 is
# not passing on three disjoint models and the zero-fill is exercised.
MLSHARED=$(awk -F'\t' 'NR>1{n=($3!=0)+($4!=0)+($5!=0); if(n>=2) c++} END{print c+0}' "$WDIR/mL.matched")
MLFILLED=$(awk -F'\t' 'NR>1{if($3==0||$4==0||$5==0) c++} END{print c+0}' "$WDIR/mL.matched")
if [ "$MLSHARED" -gt 0 ] && [ "$MLFILLED" -gt 0 ]; then
  echo "PASS  --pgs-models-list: union has shared ($MLSHARED) and zero-filled ($MLFILLED) rows"
else
  echo "FAIL  --pgs-models-list: shared=$MLSHARED zero-filled=$MLFILLED, both must exceed 0"; exit 1
fi

# Case 45: the sharp one.  Two columns carrying the same variant were built
# from the same source row, so their harmonised betas must be identical.
MLDIS=$(awk -F'\t' 'NR>1{
    if ($3!=0 && $4!=0 && $3!=$4) c++
    if ($4!=0 && $5!=0 && $4!=$5) c++
    if ($3!=0 && $5!=0 && $3!=$5) c++
  } END{print c+0}' "$WDIR/mL.matched")
if [ "$MLDIS" -eq 0 ]; then
  echo "PASS  --pgs-models-list: shared variants agree across model columns"
else
  echo "FAIL  --pgs-models-list: $MLDIS shared variants disagree between columns"; exit 1
fi

# Case 46: a one-model list reproduces --pgs-model exactly.  This is the
# invariant pinning the shared harmonisation core -- --pgs-model is the
# one-model case of the same code, so only the column header may differ.
printf '%s\n' "$WDIR/mL_A.txt" > "$WDIR/mL_one.txt"
bin/adjutant --bim "$BIM" --pgs-models-list "$WDIR/mL_one.txt" 1 2 3 4 5 \
  --out "$WDIR/mL_one" > "$WDIR/mL_one_out.txt" 2>&1
bin/adjutant --bim "$BIM" --pgs-model "$WDIR/mL_A.txt" 1 2 3 4 5 \
  --out "$WDIR/mL_single" > "$WDIR/mL_single_out.txt" 2>&1
tail -n +2 "$WDIR/mL_single.matched" > "$WDIR/mL_body_single.txt"
tail -n +2 "$WDIR/mL_one.matched"    > "$WDIR/mL_body_one.txt"
if cmp -s "$WDIR/mL_body_single.txt" "$WDIR/mL_body_one.txt"; then
  echo "PASS  --pgs-models-list: a one-model list equals --pgs-model"
else
  echo "FAIL  --pgs-models-list: one-model list differs from --pgs-model:"
  diff "$WDIR/mL_body_single.txt" "$WDIR/mL_body_one.txt" | head -5
  exit 1
fi

# Column identity.  Every assertion above is invariant to a PERMUTATION of the
# model columns -- "shared variants agree" and the row/zero counts all survive
# cyclically shifting each model's values one column along, which was verified
# by mutating the matrix write and watching the suite stay green.  This ties
# each column to ITS model's data: the mL_A column must equal, variant for
# variant, what --pgs-model produces for mL_A alone.
MLCOL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="mL_A") print i}' "$WDIR/mL.matched")
MISMATCH=$(awk -F'\t' -v c="$MLCOL" '
  NR==FNR && FNR>1 { want[$1]=$3; next }
  FNR>1 && ($1 in want) { if ($c != want[$1]) n++ }
  END { print n+0 }' "$WDIR/mL_single.matched" "$WDIR/mL.matched")
if [ -n "$MLCOL" ] && [ "$MISMATCH" -eq 0 ]; then
  echo "PASS  --pgs-models-list: each column holds its own model's betas"
else
  echo "FAIL  --pgs-models-list: mL_A column (idx '$MLCOL') differs from --pgs-model in $MISMATCH rows"
  exit 1
fi

# Case 47: the report keeps per-model provenance.
MLMRH=$(head -1 "$WDIR/mL.match_report" | cut -f1)
MLMRN=$(awk -F'\t' 'NR>1{c[$1]=1} END{print length(c)}' "$WDIR/mL.match_report")
if [ "$MLMRH" = "#MODEL" ] && [ "$MLMRN" -eq 3 ]; then
  echo "PASS  --pgs-models-list: match_report has a MODEL column naming all 3 models"
else
  echo "FAIL  --pgs-models-list: report header '$MLMRH', $MLMRN distinct models"; exit 1
fi

# Case 48: two paths yielding one name is refused, not silently disambiguated
# -- the names become column headers the user addresses with --score-col-nums.
mkdir -p "$WDIR/mLsub"
cp "$WDIR/mL_A.txt" "$WDIR/mLsub/mL_A.txt"
printf '%s\n%s\n' "$WDIR/mL_A.txt" "$WDIR/mLsub/mL_A.txt" > "$WDIR/mL_dup.txt"
if bin/adjutant --bim "$BIM" --pgs-models-list "$WDIR/mL_dup.txt" 1 2 3 4 5 \
     --out "$WDIR/mL_dup" > "$WDIR/mL_dup_out.txt" 2>&1
then
  echo "FAIL  --pgs-models-list: duplicate model names were accepted"; exit 1
else
  if saidflat "$WDIR/mL_dup_out.txt" "both give the model name"; then
    echo "PASS  --pgs-models-list: duplicate model names refused"
  else
    echo "FAIL  --pgs-models-list: duplicates refused with the wrong message"
    cat "$WDIR/mL_dup_out.txt"; exit 1
  fi
fi

# Case 49: a model matching nothing warns and loses its column, but the run
# survives, because the other models are still usable.
awk 'BEGIN{OFS="\t"} {$2=$2+90000000; print}' "$WDIR/mL_A.txt" > "$WDIR/mL_none.txt"
printf '%s\n%s\n' "$WDIR/mL_A.txt" "$WDIR/mL_none.txt" > "$WDIR/mL_nl.txt"
bin/adjutant --bim "$BIM" --pgs-models-list "$WDIR/mL_nl.txt" 1 2 3 4 5 \
  --out "$WDIR/mL_nm" > "$WDIR/mL_nm_out.txt" 2>&1
MLNMH=$(head -1 "$WDIR/mL_nm.matched")
if saidflat "$WDIR/mL_nm_out.txt" "matched no variants" &&
   [ "$MLNMH" = "$(printf 'ID\tA1\tmL_A')" ]; then
  echo "PASS  --pgs-models-list: a model matching nothing warns and drops its column"
else
  echo "FAIL  --pgs-models-list: header '$MLNMH' after a zero-match model"
  cat "$WDIR/mL_nm_out.txt"; exit 1
fi

# Cases 50-52: the guards.  Each must be refused, with a message naming why.
ml_refuse() {  # ml_refuse <label> <expected substring> <args...>
  _lbl=$1; _want=$2; shift 2
  if bin/adjutant "$@" --out "$WDIR/mL_rej" > "$WDIR/mL_rej_out.txt" 2>&1; then
    echo "FAIL  --pgs-models-list: $_lbl was accepted"; exit 1
  fi
  if saidflat "$WDIR/mL_rej_out.txt" "$_want"; then
    echo "PASS  --pgs-models-list: refused $_lbl"
  else
    echo "FAIL  --pgs-models-list: $_lbl refused without '$_want'"
    cat "$WDIR/mL_rej_out.txt"; exit 1
  fi
}
ml_refuse "without --bim" "requires --bim" \
  --pgs-models-list "$WDIR/mL_list.txt" 1 2 3 4 5
ml_refuse "alongside --pgs-model" "cannot be used with --pgs-models-list" \
  --bim "$BIM" --pgs-models-list "$WDIR/mL_list.txt" 1 2 3 4 5 \
  --pgs-model "$WDIR/mL_A.txt" 1 2 3 4 5
ml_refuse "a list with no paths" "no model paths" \
  --bim "$BIM" --pgs-models-list "$WDIR/mL_empty.txt" 1 2 3 4 5

# ---------------------------------------------------------------------------
# Cases 58-59: harmonisation progress output.
#
# Matching is the phase with no output of its own, and on a large weight file
# against a large target it is where the wall-clock goes.  It now reports
# progress the same way "Calculating allele frequencies..." does: the text via
# logputs so the log keeps it, the percentage written straight to stdout and
# backspaced away.
#
# The second case is the one that matters.  Emitting the percentage through
# logprintf would look identical on a terminal -- the backspaces still erase it
# on screen -- while quietly filling the .log with 0%1%2%...99% on every run.
# Only inspecting the log catches that.
# ---------------------------------------------------------------------------
bin/adjutant --bim "$BIM" --pgs-model "$WDIR/pgs_src.txt" 1 2 3 4 5 \
  --out "$WDIR/prog" > "$WDIR/prog_out.txt" 2>&1 || true

# Deliberately tolerant of what sits between: if the percentage counter ever
# leaks into the log it lands right here, and that is the NEXT case's job to
# catch.  A stricter pattern would make this case fail first and mask it.
if saidflat "$WDIR/prog.log" "matching .*done\."; then
  echo "PASS  --pgs-model: reports matching progress, and completion"
else
  echo "FAIL  --pgs-model: no 'matching ... done.' line in the log"
  grep -c . "$WDIR/prog.log"; exit 1
fi

PCT_IN_LOG=$(grep -c "%" "$WDIR/prog.log" || true)
if [ "$PCT_IN_LOG" -eq 0 ]; then
  echo "PASS  --pgs-model: the percentage counter stays out of the .log"
else
  echo "FAIL  --pgs-model: $PCT_IN_LOG log lines contain '%'; the counter is being logged"
  grep -m2 "%" "$WDIR/prog.log"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 60: the old --resolve-ambig spelling is gone.
#
# The suite exercises --resolve-ambiguous throughout, so the new spelling is
# well covered by construction.  Nothing covered the other half, and that is
# the half that broke during the rename: the parser matches on flagname_p2 --
# the name with its first character already consumed by the switch -- so the
# token in the source reads "esolve-ambig" and a sweep for "resolve-ambig"
# leaves it untouched.  Everything then builds, every message and the help say
# --resolve-ambiguous, and the parser still accepts only the old spelling.
# ---------------------------------------------------------------------------
if bin/adjutant --bim "$BIM" --pgs-model "$WDIR/pgs_src.txt" 1 2 3 4 5 \
     --resolve-ambig --out "$WDIR/oldspell" > "$WDIR/oldspell_cli.txt" 2>&1
then
  echo "FAIL  --resolve-ambig: old spelling still accepted"; exit 1
elif saidflat "$WDIR/oldspell_cli.txt" "Unrecognized flag ('--resolve-ambig')"; then
  echo "PASS  --resolve-ambiguous: the old --resolve-ambig spelling is rejected"
else
  echo "FAIL  --resolve-ambig: refused, but not as an unrecognized flag"
  cat "$WDIR/oldspell_cli.txt"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 61-69: --give-corr-matrix.
#
# The numbers themselves were adjudicated once against R:
#
#   d <- read.table("x.sscore", header=TRUE, comment.char="", sep="\t")
#   r <- cor(as.matrix(d[, grep("_AVG$", names(d))]))
#   a <- as.matrix(read.table("x.pgs_corr", header=TRUE, comment.char="",
#                             sep="\t", row.names=1))
#   max(abs(r - a))                 # 2.3e-13
#
# That run is not repeated here: the suite must stay runnable without R, and a
# Pearson correlation is short enough to recompute honestly in awk.  What the
# awk reference below is actually for is the class of bug that does occur --
# wrong column pairing, a transposed matrix, the wrong sample set, the wrong
# denominator -- not the formula itself, which R has already settled.
# ---------------------------------------------------------------------------

# cell <file> <ROW NAME> <COL NAME> -> one matrix entry, addressed by name in
# BOTH dimensions.
#
# Positional indexing is worse than useless on this file: a correlation matrix
# is symmetric with a constant diagonal, so a transposed matrix, or one whose
# columns were permuted consistently, is indistinguishable from the right one
# under $2/$3.  Only the names distinguish them.
cell() {
  _v=$(awk -F'\t' -v r="$2" -v c="$3" '
    NR==1 { for (i = 2; i <= NF; i++) if ($i == c) ci = i; next }
    ci && $1 == r { print $ci; exit }' "$1")
  [ -n "$_v" ] || { echo "FAIL  $1 has no cell [$2, $3]" >&2; exit 1; }
  printf '%s\n' "$_v"
}

# pearson <sscore file> <COL A> <COL B> -> Pearson r straight from the .sscore.
#
# Subscripts run 1..m.  Written 0-based first -- `x[m] = ...; m++` with m
# unset -- and awk converts an uninitialised subscript to "" rather than 0, so
# row 1 landed in x[""] and the END loop read an unset x[0] as a real 0.  The
# reference then correlated 2503 samples against one fabricated zero and came
# out wrong in the 4th decimal: large enough to fail, small enough to look like
# a genuine disagreement with the C rather than a bug in the check.
pearson() {
  awk -F'\t' -v a="$2" -v b="$3" '
    NR==1 { for (i = 1; i <= NF; i++) { n = $i; sub(/^#/, "", n)
              if (n == a) ia = i; if (n == b) ib = i }
            if (!ia || !ib) exit 1
            next }
    { m++; x[m] = $ia + 0; y[m] = $ib + 0 }
    END { if (m < 2) exit 1
          for (i = 1; i <= m; i++) { sx += x[i]; sy += y[i] }
          mx = sx / m; my = sy / m
          for (i = 1; i <= m; i++) { dx = x[i] - mx; dy = y[i] - my
            sxy += dx * dy; sxx += dx * dx; syy += dy * dy }
          printf "%.12g\n", sxy / sqrt(sxx * syy) }' "$1"
}

# close <a> <b> <tol> -> exit 0 if |a-b| <= tol
close() {
  awk -v a="$1" -v b="$2" -v t="$3" 'BEGIN { exit !((a - b) ** 2 <= t * t) }'
}



# `|| true` on every run below, and it is load-bearing, not decoration.  The
# script runs under `set -e`: a build that refuses one of these commands would
# abort the whole suite here, silently, before a single assertion ran -- and a
# red-check of the guards in Case 69 reported "no failure" for exactly that
# reason while the mutation was working perfectly.  Let the command fail and
# let the assertion say so.
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --give-corr-matrix --out "$WDIR/cm" > "$WDIR/cm_cli.txt" 2>&1 || true
need "$WDIR/cm.pgs_corr" "--give-corr-matrix: the base --sscore run produced no matrix"

# ---------------------------------------------------------------------------
# Case 61: the header is exactly this.  The one place .pgs_corr's layout is
# stated rather than looked up -- cell() resolves everything else by name, so
# without this the names themselves are unchecked.
# ---------------------------------------------------------------------------
CM_EXPECTED_HEADER='#ID	W1	W2	W3'
CM_ACTUAL_HEADER=$(head -1 "$WDIR/cm.pgs_corr")
if [ "$CM_ACTUAL_HEADER" = "$CM_EXPECTED_HEADER" ]; then
  echo "PASS  --give-corr-matrix: .pgs_corr header exact (names and order)"
else
  echo "FAIL  --give-corr-matrix: .pgs_corr header mismatch"
  echo "  expected: $CM_EXPECTED_HEADER"
  echo "  got:      $CM_ACTUAL_HEADER"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 62: full square.  N models give N+1 lines of N+1 fields each -- the
# property that lets a reader skip knowing which triangle it got.
# ---------------------------------------------------------------------------
CM_N=$(awk -F'\t' 'NR==1{print NF-1; exit}' "$WDIR/cm.pgs_corr")
CM_BAD=$(awk -F'\t' -v n="$CM_N" '
  NF != n + 1 { print NR": "NF" fields" }
  END { if (NR != n + 1) print "line count "NR", expected "n+1 }' "$WDIR/cm.pgs_corr")
if [ "$CM_N" -eq 3 ] && [ -z "$CM_BAD" ]; then
  echo "PASS  --give-corr-matrix: 3 models give a full 3x3 square"
else
  echo "FAIL  --give-corr-matrix: not a full square ($CM_N models):"
  printf '%s\n' "$CM_BAD"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 63: every off-diagonal entry matches the correlation computed from the
# .sscore independently, to 1e-10.
#
# Addressed by name, and the two directions are read separately rather than one
# being assumed from the other -- symmetry is Case 64's claim, not this one's.
# ---------------------------------------------------------------------------
CM_BADVAL=""
for PAIR in "W1 W2" "W1 W3" "W2 W3"; do
  set -- $PAIR
  EXPECT=$(pearson "$WDIR/direct.sscore" "$1_AVG" "$2_AVG")
  [ -n "$EXPECT" ] || { echo "FAIL  --give-corr-matrix: no reference for $1/$2"; exit 1; }
  for CELL in "$1 $2" "$2 $1"; do
    set -- $CELL
    GOT=$(cell "$WDIR/cm.pgs_corr" "$1" "$2")
    close "$GOT" "$EXPECT" 1e-10 || CM_BADVAL="$CM_BADVAL [$1,$2] got $GOT want $EXPECT"
  done
done
if [ -z "$CM_BADVAL" ]; then
  echo "PASS  --give-corr-matrix: off-diagonal entries match an independent Pearson to 1e-10"
else
  echo "FAIL  --give-corr-matrix:$CM_BADVAL"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 64: the diagonal is the literal string "1".
#
# The pair loop in the writer starts one past the diagonal, so nothing computes
# these entries -- they are assigned.  That makes them exactly the kind of value
# that can be left at whatever the buffer held: verified by mutation, assigning
# 0 instead fails here on all three.
#
# A literal string comparison rather than a tolerance, because there is no
# arithmetic to tolerate.  (Computing the diagonal through the general formula
# was tried too, and prints "1" as well -- so this case is pinning that the
# entry is written at all, not defending a numerical short-circuit.)
# ---------------------------------------------------------------------------
CM_BADDIAG=""
for M in W1 W2 W3; do
  D=$(cell "$WDIR/cm.pgs_corr" "$M" "$M")
  [ "$D" = "1" ] || CM_BADDIAG="$CM_BADDIAG [$M,$M]=$D"
done
if [ -z "$CM_BADDIAG" ]; then
  echo "PASS  --give-corr-matrix: the diagonal is exactly 1"
else
  echo "FAIL  --give-corr-matrix: diagonal is not exactly 1:$CM_BADDIAG"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 65: the matrix equals its own transpose, byte for byte.
# ---------------------------------------------------------------------------
awk -F'\t' '{ for (i = 1; i <= NF; i++) t[i, NR] = $i; nr = NR; nf = NF }
  END { for (i = 1; i <= nf; i++) { s = t[i, 1]
          for (j = 2; j <= nr; j++) s = s "\t" t[i, j]
          print s } }' "$WDIR/cm.pgs_corr" > "$WDIR/cm_transposed.txt"
if diff "$WDIR/cm.pgs_corr" "$WDIR/cm_transposed.txt" > /dev/null 2>&1; then
  echo "PASS  --give-corr-matrix: the matrix is symmetric"
else
  echo "FAIL  --give-corr-matrix: matrix differs from its transpose:"
  diff "$WDIR/cm.pgs_corr" "$WDIR/cm_transposed.txt" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Case 66: "a pure property of the scores", stated as something that can fail.
#
# The same scores with a phenotype and covariates attached must give a
# byte-identical matrix -- and the phenotype used here has missing values, so
# the complete-case set --give-best-pgs fits on is genuinely smaller than the
# score matrix.  Without the missingness the two sample sets would coincide and
# this case would pass whether or not the restriction was being applied.
# ---------------------------------------------------------------------------
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/missing.txt" --pheno-name MISSQT \
  --covar "$WDIR/covar.txt" \
  --give-best-pgs --give-corr-matrix --out "$WDIR/cm_pheno" > "$WDIR/cm_pheno_cli.txt" 2>&1 || true
need "$WDIR/cm_pheno.pgs_corr" "--give-corr-matrix: the run with --pheno/--covar produced no matrix"

CM_FIT_N=$(grep -v '^#' "$WDIR/cm_pheno.best_pgs" | awk -v c="$C_N" '{print $c}' | sort -u)
CM_ALL_N=$(($(wc -l < "$WDIR/direct.sscore") - 1))
if [ "$CM_FIT_N" = "$CM_ALL_N" ]; then
  echo "FAIL  --give-corr-matrix: the fitted set ($CM_FIT_N) equals the score set;"
  echo "      case 66 cannot distinguish the two sample sets"; exit 1
elif diff "$WDIR/cm.pgs_corr" "$WDIR/cm_pheno.pgs_corr" > /dev/null 2>&1; then
  echo "PASS  --give-corr-matrix: identical with and without --pheno/--covar (fit on $CM_FIT_N of $CM_ALL_N)"
else
  echo "FAIL  --give-corr-matrix: --pheno/--covar changed the matrix:"
  diff "$WDIR/cm.pgs_corr" "$WDIR/cm_pheno.pgs_corr" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Case 67: both routes, same scores, same matrix.
#
# --pfile computes the scores and correlates them in one run; --sscore reads
# that same .sscore back.  One shared helper is meant to make these identical,
# so anything that reaches only one of the two call sites shows up here.
# ---------------------------------------------------------------------------
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-corr-matrix --out "$WDIR/cm_pfile" > "$WDIR/cm_pfile_cli.txt" 2>&1 || true
need "$WDIR/cm_pfile.pgs_corr" "--give-corr-matrix: the --pfile route produced no matrix"

if diff "$WDIR/cm.pgs_corr" "$WDIR/cm_pfile.pgs_corr" > /dev/null 2>&1; then
  echo "PASS  --give-corr-matrix: --pfile and --sscore routes agree byte for byte"
else
  echo "FAIL  --give-corr-matrix: the two routes disagree:"
  diff "$WDIR/cm.pgs_corr" "$WDIR/cm_pfile.pgs_corr" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Case 68: a zero-variance column is NA across its whole row and column.
#
# Not 0.  0 reads as "uncorrelated", which is a claim about two things that
# vary; a constant column has no correlation with anything, and the file has to
# say which of those two it means.  The surviving 3x3 block must be untouched,
# so a degenerate column cannot quietly poison its neighbours.
# ---------------------------------------------------------------------------
awk -F'\t' 'NR==1{print $0"\tWCONST_AVG"} NR>1{print $0"\t5"}' \
  "$WDIR/direct.sscore" > "$WDIR/const.sscore"
bin/adjutant --sscore "$WDIR/const.sscore" \
  --give-corr-matrix --out "$WDIR/cm_const" > "$WDIR/cm_const_cli.txt" 2>&1 || true
need "$WDIR/cm_const.pgs_corr" "--give-corr-matrix: the constant-column run produced no matrix"

CM_CONST_BAD=""
for M in W1 W2 W3 WCONST; do
  for DIR in "WCONST $M" "$M WCONST"; do
    set -- $DIR
    V=$(cell "$WDIR/cm_const.pgs_corr" "$1" "$2")
    [ "$V" = "NA" ] || CM_CONST_BAD="$CM_CONST_BAD [$1,$2]=$V"
  done
done
# Strike out WCONST's row and its column, by name, and what is left must be the
# matrix from before it was added.
awk -F'\t' -v drop=WCONST '
  NR == 1 { for (i = 2; i <= NF; i++) if ($i == drop) di = i }
  $1 == drop { next }
  { s = $1; for (i = 2; i <= NF; i++) if (i != di) s = s "\t" $i; print s }' \
  "$WDIR/cm_const.pgs_corr" > "$WDIR/cm_const_block.txt"

if [ -n "$CM_CONST_BAD" ]; then
  echo "FAIL  --give-corr-matrix: constant column is not NA:$CM_CONST_BAD"; exit 1
elif ! saidflat "$WDIR/cm_const_cli.txt" "1 score column is constant or non-finite"; then
  echo "FAIL  --give-corr-matrix: constant column produced no warning:"
  cat "$WDIR/cm_const_cli.txt"; exit 1
elif ! diff "$WDIR/cm_const_block.txt" "$WDIR/cm.pgs_corr" > /dev/null 2>&1; then
  echo "FAIL  --give-corr-matrix: a constant column changed the other entries:"
  diff "$WDIR/cm_const_block.txt" "$WDIR/cm.pgs_corr" | head -6; exit 1
else
  echo "PASS  --give-corr-matrix: a constant column is NA, warned about, and inert"
fi

# ---------------------------------------------------------------------------
# Case 69: the guards.
#
# The relaxation is the point: --sscore no longer requires --give-best-pgs, and
# --give-corr-matrix alone no longer requires --pheno.  Both used to be refused,
# so a regression reads as a hard error rather than a wrong number.
#
# The first line duplicates a precondition the base run at the top of this block
# already enforces -- any mutation that refuses it takes out the need() there
# first.  It is kept as the explicit statement of the relaxation, but the base
# run is what actually guards it; the other three lines fail on their own terms
# and were each seen to.
# ---------------------------------------------------------------------------
cm_expect() {  # cm_expect <label> <"ok"|message fragment> <args...>
  _label=$1; _want=$2; shift 2
  if bin/adjutant "$@" > "$WDIR/cm_g.txt" 2>&1; then
    if [ "$_want" = "ok" ]; then
      echo "PASS  --give-corr-matrix: $_label"
    else
      echo "FAIL  --give-corr-matrix: $_label -- accepted, expected refusal"; exit 1
    fi
  elif [ "$_want" = "ok" ]; then
    echo "FAIL  --give-corr-matrix: $_label -- refused, expected success:"
    cat "$WDIR/cm_g.txt"; exit 1
  elif saidflat "$WDIR/cm_g.txt" "$_want"; then
    echo "PASS  --give-corr-matrix: $_label"
  else
    echo "FAIL  --give-corr-matrix: $_label -- refused, but not with '$_want':"
    cat "$WDIR/cm_g.txt"; exit 1
  fi
}
cm_expect "--sscore without --pheno is now accepted" ok \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --give-corr-matrix --out "$WDIR/cm_g1"
cm_expect "--sscore with neither report is refused" \
  "requires --give-best-pgs or --give-corr-matrix" \
  --sscore "$WDIR/direct.sscore" --out "$WDIR/cm_g2"
cm_expect "without --score or --sscore it is refused, naming itself" \
  "give-corr-matrix must be used with --score or --sscore" \
  --pfile "$PFILE_PFX" vzs --give-corr-matrix --out "$WDIR/cm_g3"
cm_expect "--sscore with --give-best-pgs still requires --pheno" \
  "requires --pheno" \
  --sscore "$WDIR/direct.sscore" --give-best-pgs --out "$WDIR/cm_g4"

# ---------------------------------------------------------------------------
# Cases 70-79: --pgs-pheno-pairs.
#
# Three phenotypes and three score columns, so the unrestricted run is nine
# rows and a two-pair file has eight ways to be wrong while still producing two
# rows.  Every case below therefore checks row IDENTITY -- the (PGS_ID, PHENO)
# set -- and never the row count.  "It returned two rows" is not the claim.
#
# `|| true` plus need() on every run, for the reason spelled out in the
# correlation block above: under set -e a refused run takes the suite out
# before any assertion prints.
# ---------------------------------------------------------------------------
awk 'NR==1{print "#IID\tPQT1\tPQT2\tPCC1"}
     NR>1{print $1"\t"(NR%97)/13.0"\t"(NR%53)/7.0"\t"((NR%3==0)?2:1)}' \
  "$PSAM" > "$WDIR/pp_pheno.txt"

# The (PGS_ID, PHENO) set of a .best_pgs, sorted and flattened to one line.
pairset() {
  awk -F'\t' -v a="$C_PGS_ID" -v b="$C_PHENO" 'NR>1{print $a"/"$b}' "$1" \
    | sort | tr '\n' ' '
}

bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --covar "$WDIR/covar.txt" \
  --give-best-pgs --out "$WDIR/pp_all" > "$WDIR/pp_all_cli.txt" 2>&1 || true
need "$WDIR/pp_all.best_pgs" "--pgs-pheno-pairs: the unrestricted run produced nothing"

PP_ALL_SET=$(pairset "$WDIR/pp_all.best_pgs")
PP_EXPECT_ALL='W1/PCC1 W1/PQT1 W1/PQT2 W2/PCC1 W2/PQT1 W2/PQT2 W3/PCC1 W3/PQT1 W3/PQT2 '
if [ "$PP_ALL_SET" = "$PP_EXPECT_ALL" ]; then
  echo "PASS  --pgs-pheno-pairs: the unrestricted run is the full 3x3 cross-product"
else
  echo "FAIL  --pgs-pheno-pairs: unrestricted run is not the cross-product"
  echo "  expected: $PP_EXPECT_ALL"
  echo "  got:      $PP_ALL_SET"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 70: the right rows, not merely the right number of them.
#
# Comments and a trailing blank line are in the file on purpose -- they are
# part of the format, and a parser that mishandled either would change the
# pair set rather than fail loudly.
# ---------------------------------------------------------------------------
printf '# phenotype  PGS\nPQT1\tW1\nPCC1\tW3\nPQT1\tW2\n\n' > "$WDIR/pp_three.txt"
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --covar "$WDIR/covar.txt" \
  --give-best-pgs --pgs-pheno-pairs "$WDIR/pp_three.txt" \
  --out "$WDIR/pp_sub" > "$WDIR/pp_sub_cli.txt" 2>&1 || true
need "$WDIR/pp_sub.best_pgs" "--pgs-pheno-pairs: the filtered run produced nothing"

PP_SUB_SET=$(pairset "$WDIR/pp_sub.best_pgs")
PP_EXPECT_SUB='W1/PQT1 W2/PQT1 W3/PCC1 '
if [ "$PP_SUB_SET" = "$PP_EXPECT_SUB" ]; then
  echo "PASS  --pgs-pheno-pairs: exactly the listed pairs, by identity"
else
  echo "FAIL  --pgs-pheno-pairs: wrong pairs"
  echo "  expected: $PP_EXPECT_SUB"
  echo "  got:      $PP_SUB_SET"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 71: the invariant that pins the filter as a pure restriction.
#
# A pairs file naming the whole cross-product must reproduce the unrestricted
# .best_pgs byte for byte -- same rows, same order, same digits.  Anything that
# made the filtered path a second implementation rather than a subset of the
# first shows up here and nowhere else.
# ---------------------------------------------------------------------------
: > "$WDIR/pp_cross.txt"
for PPH in PQT1 PQT2 PCC1; do
  for PSC in W1 W2 W3; do
    printf '%s\t%s\n' "$PPH" "$PSC" >> "$WDIR/pp_cross.txt"
  done
done
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --covar "$WDIR/covar.txt" \
  --give-best-pgs --pgs-pheno-pairs "$WDIR/pp_cross.txt" \
  --out "$WDIR/pp_full" > "$WDIR/pp_full_cli.txt" 2>&1 || true
need "$WDIR/pp_full.best_pgs" "--pgs-pheno-pairs: the cross-product run produced nothing"

if diff "$WDIR/pp_all.best_pgs" "$WDIR/pp_full.best_pgs" > /dev/null 2>&1; then
  echo "PASS  --pgs-pheno-pairs: the full cross-product reproduces the unrestricted output byte for byte"
else
  echo "FAIL  --pgs-pheno-pairs: naming every pair changed the output:"
  diff "$WDIR/pp_all.best_pgs" "$WDIR/pp_full.best_pgs" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Case 72: the output does not depend on how the pairs file was written.
#
# Same three pairs reversed, with a duplicate added.  Rows are emitted in
# score-column order and duplicates collapse, so the file must be identical --
# not merely the same set.
# ---------------------------------------------------------------------------
printf 'PQT1\tW2\nPCC1\tW3\nPQT1\tW1\nPCC1\tW3\n' > "$WDIR/pp_perm.txt"
bin/adjutant \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --covar "$WDIR/covar.txt" \
  --give-best-pgs --pgs-pheno-pairs "$WDIR/pp_perm.txt" \
  --out "$WDIR/pp_perm" > "$WDIR/pp_perm_cli.txt" 2>&1 || true
need "$WDIR/pp_perm.best_pgs" "--pgs-pheno-pairs: the reordered run produced nothing"

if diff "$WDIR/pp_sub.best_pgs" "$WDIR/pp_perm.best_pgs" > /dev/null 2>&1; then
  echo "PASS  --pgs-pheno-pairs: order and duplicates in the pairs file do not change the output"
else
  echo "FAIL  --pgs-pheno-pairs: reordering or duplicating pairs changed the output:"
  diff "$WDIR/pp_sub.best_pgs" "$WDIR/pp_perm.best_pgs" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Case 73: a phenotype paired with several models is still ranked.
#
# The filter removes rows; it must not disturb the sort within what is left.
# PQT1 keeps two of the three models, and its block must still be R2_INCR
# descending.
# ---------------------------------------------------------------------------
PP_ORDER_BAD=$(awk -F'\t' -v p="$C_PHENO" -v r="$C_R2_INCR" '
  NR > 1 && $p == "PQT1" { if (seen && $r > prev + 1e-12) print "row "NR" out of order"
                           prev = $r; seen = 1 }
  END { if (!seen) print "no PQT1 rows at all" }' "$WDIR/pp_sub.best_pgs")
PP_QT1_ROWS=$(awk -F'\t' -v p="$C_PHENO" 'NR>1 && $p=="PQT1"' "$WDIR/pp_sub.best_pgs" | wc -l | tr -d ' ')
if [ -z "$PP_ORDER_BAD" ] && [ "$PP_QT1_ROWS" = "2" ]; then
  echo "PASS  --pgs-pheno-pairs: a multi-pair phenotype keeps its R2_INCR ordering"
else
  echo "FAIL  --pgs-pheno-pairs: PQT1 has $PP_QT1_ROWS rows; $PP_ORDER_BAD"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 74: both routes, same pairs, same rows.
#
# The filter lives at FitWritePhenoBlock precisely so there is one of it.  This
# is what would catch a second copy appearing on one route.
# ---------------------------------------------------------------------------
bin/adjutant \
  --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --pheno "$WDIR/pp_pheno.txt" --covar "$WDIR/covar.txt" \
  --give-best-pgs --pgs-pheno-pairs "$WDIR/pp_three.txt" \
  --out "$WDIR/pp_pfile" > "$WDIR/pp_pfile_cli.txt" 2>&1 || true
need "$WDIR/pp_pfile.best_pgs" "--pgs-pheno-pairs: the --pfile run produced nothing"

grep -v '^#' "$WDIR/pp_sub.best_pgs"   > "$WDIR/pp_sub_data.txt"
grep -v '^#' "$WDIR/pp_pfile.best_pgs" > "$WDIR/pp_pfile_data.txt"
if diff "$WDIR/pp_sub_data.txt" "$WDIR/pp_pfile_data.txt" > /dev/null 2>&1; then
  echo "PASS  --pgs-pheno-pairs: --pfile and --sscore routes filter identically"
else
  echo "FAIL  --pgs-pheno-pairs: the two routes disagree under a filter:"
  diff "$WDIR/pp_sub_data.txt" "$WDIR/pp_pfile_data.txt" | head -6; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 75-79: the guards.  A missing name is fatal rather than skipped --
# a typo that silently dropped a pair would leave a report that looks complete
# and is not, which is the failure this flag exists to prevent.
# ---------------------------------------------------------------------------
pp_refuse() {  # pp_refuse <label> <message fragment> <args...>
  _label=$1; _want=$2; shift 2
  if bin/adjutant "$@" > "$WDIR/pp_rej.txt" 2>&1; then
    echo "FAIL  --pgs-pheno-pairs: $_label -- accepted, expected refusal"; exit 1
  elif saidflat "$WDIR/pp_rej.txt" "$_want"; then
    echo "PASS  --pgs-pheno-pairs: refused $_label"
  else
    echo "FAIL  --pgs-pheno-pairs: $_label -- refused, but not with '$_want':"
    cat "$WDIR/pp_rej.txt"; exit 1
  fi
}
printf 'NOSUCHPHENO\tW1\n'  > "$WDIR/pp_badpheno.txt"
printf 'PQT1\tNOSUCHMODEL\n' > "$WDIR/pp_badpgs.txt"
printf 'PQT1\n'              > "$WDIR/pp_onecol.txt"
printf '# only a comment\n'  > "$WDIR/pp_empty.txt"
printf 'PQT1\tW1\n\nPCC1\tW3\n' > "$WDIR/pp_blank.txt"

pp_refuse "an unknown phenotype, naming it" "'NOSUCHPHENO' is not a phenotype" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_badpheno.txt" --out "$WDIR/pp_g1"
pp_refuse "an unknown PGS, naming it" "'NOSUCHMODEL' is not a score column" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_badpgs.txt" --out "$WDIR/pp_g2"
pp_refuse "a one-column line, naming the line number" "line 1 of .* has one column" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_onecol.txt" --out "$WDIR/pp_g3"
pp_refuse "a file with no pairs" "contains no pairs" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_empty.txt" --out "$WDIR/pp_g4"
pp_refuse "--pheno-name alongside it, naming both flags" \
  "cannot be used with --pheno-name" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --pheno-name PQT1 --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_three.txt" --out "$WDIR/pp_g5"
pp_refuse "--pheno-col-nums alongside it, naming that flag instead" \
  "cannot be used with --pheno-col-nums" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --pheno-col-nums 1 --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_three.txt" --out "$WDIR/pp_g6"
pp_refuse "use without --give-best-pgs" "must be used with --give-best-pgs" \
  --pfile "$PFILE_PFX" vzs --score "$WDIR/w3.txt" 1 2 header-read \
  --pgs-pheno-pairs "$WDIR/pp_three.txt" --out "$WDIR/pp_g7"

# An interior blank line gets adjutant's own message, not plink2's internal
# "Unexpected interior empty line" -- TextStream stops AT the blank and only
# then reports content past it, so without this the pair list would truncate.
pp_refuse "an interior blank line, in its own words" "blank line inside" \
  --sscore "$WDIR/direct.sscore" --score-col-nums "$SCORE_COL_RANGE" \
  --pheno "$WDIR/pp_pheno.txt" --give-best-pgs \
  --pgs-pheno-pairs "$WDIR/pp_blank.txt" --out "$WDIR/pp_g8"

# ---------------------------------------------------------------------------
# Cases 80-89: --score-partitions.
#
# Worth stating plainly: this flag is the FIRST live user of the
# --q-score-range machinery inside ScoreReport.  adjutant has no
# --q-score-range parser case, so qsr_range_fname is permanently null and that
# whole path -- per-partition denominators in allele_ct_bases[], the
# [sample x partition x model] accumulation, the per-partition missingness
# accounting -- has never executed in any test in this repository.  It is
# effectively new code.  The --extract equivalence below is therefore not a
# formality; it is the first thing that has ever checked those denominators.
# ---------------------------------------------------------------------------
# PHALF and POTHER tile the variant set exactly (needed for additivity);
# PSMALL is a sparse overlapping subset.
awk 'NR==1{print "SNP\tPHALF\tPOTHER\tPSMALL"; next}
     {print $1"\t"(NR%2)"\t"(1-NR%2)"\t"((NR%5==0)?1:0)}' "$WDIR/w3.txt" > "$WDIR/sp_parts.txt"

bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read cols=+scoresums --score-col-nums 3-5 \
  --score-partitions "$WDIR/sp_parts.txt" \
  --out "$WDIR/sp" > "$WDIR/sp_cli.txt" 2>&1 || true
need "$WDIR/sp.sscore" "--score-partitions: the partitioned run produced nothing"

# ---------------------------------------------------------------------------
# Case 80: the score columns are exactly these, in this order.
#
# Pins the naming rule -- total first with no suffix, then partitions in file
# order, grouped by model -- in the one place it is stated rather than derived.
# Only the score columns are pinned: the leading block carries .psam phenotypes
# that have nothing to do with this feature.
# ---------------------------------------------------------------------------
# NAMED_ALLELE_DOSAGE_SUM also ends in _SUM and is a bookkeeping column, not a
# score -- the same collision --sscore's auto-detection has to exclude by name.
SP_COLS=$(head -1 "$WDIR/sp.sscore" | tr '\t' '\n' | grep -E '_(AVG|SUM)$' \
          | grep -v '^NAMED_ALLELE_DOSAGE_SUM$' | tr '\n' ' ')
SP_WANT='W1_AVG W1_PHALF_AVG W1_POTHER_AVG W1_PSMALL_AVG W2_AVG W2_PHALF_AVG W2_POTHER_AVG W2_PSMALL_AVG W3_AVG W3_PHALF_AVG W3_POTHER_AVG W3_PSMALL_AVG W1_SUM W1_PHALF_SUM W1_POTHER_SUM W1_PSMALL_SUM W2_SUM W2_PHALF_SUM W2_POTHER_SUM W2_PSMALL_SUM W3_SUM W3_PHALF_SUM W3_POTHER_SUM W3_PSMALL_SUM '
if [ "$SP_COLS" = "$SP_WANT" ]; then
  echo "PASS  --score-partitions: score columns exact (total first, then partitions, grouped by model)"
else
  echo "FAIL  --score-partitions: score column layout wrong"
  echo "  expected: $SP_WANT"
  echo "  got:      $SP_COLS"; exit 1
fi

# spcol <file> <COLUMN NAME> -> that column's values, one per line
spcol() {
  awk -F'\t' -v want="$2" 'NR==1{for(i=1;i<=NF;i++){n=$i; sub(/^#/,"",n); if(n==want) c=i}
                                 if(!c){exit 1}; next} {print $c}' "$1"
}

# spsame <fileA> <COLA> <fileB> <COLB> -> true if the two columns are identical.
#
# Written with temp files, not `diff <(...) <(...)`: process substitution is a
# bashism and this script is #!/bin/sh, where it is a syntax error that takes
# the whole run out silently.
#
# The emptiness check is the load-bearing part.  A missing column makes spcol
# produce nothing, and two empty outputs compare EQUAL -- so without it, a
# renamed or absent column would read as a pass rather than a failure. That is
# the vacuous-assertion shape this suite exists to avoid.
spsame() {
  spcol "$1" "$2" > "$WDIR/sp_lhs.txt" 2>/dev/null
  spcol "$3" "$4" > "$WDIR/sp_rhs.txt" 2>/dev/null
  [ -s "$WDIR/sp_lhs.txt" ] && [ -s "$WDIR/sp_rhs.txt" ] || return 2
  diff "$WDIR/sp_lhs.txt" "$WDIR/sp_rhs.txt" > /dev/null 2>&1
}

# ---------------------------------------------------------------------------
# Case 81: the total column is plain --score, to the byte.
#
# The invariant that pins partitioning as a decomposition of one scoring pass
# rather than a second implementation of it.
# ---------------------------------------------------------------------------
bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --out "$WDIR/sp_plain" > "$WDIR/sp_plain_cli.txt" 2>&1 || true
need "$WDIR/sp_plain.sscore" "--score-partitions: the unpartitioned run produced nothing"

SP_TOTBAD=""
for M in W1 W2 W3; do
  spsame "$WDIR/sp_plain.sscore" "${M}_AVG" "$WDIR/sp.sscore" "${M}_AVG" \
    || SP_TOTBAD="$SP_TOTBAD $M"
done
if [ -z "$SP_TOTBAD" ]; then
  echo "PASS  --score-partitions: every total column equals plain --score byte for byte"
else
  echo "FAIL  --score-partitions: total column differs from plain --score for:$SP_TOTBAD"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 82: each partition equals --extract restricted to its own variants.
#
# This is what makes "per-partition denominator" a checkable claim instead of a
# description: --extract genuinely rescores a smaller variant set, so its
# denominator is that set's, and the two must agree exactly.
# ---------------------------------------------------------------------------
SP_PARTBAD=""
for P in PHALF POTHER PSMALL; do
  awk -F'\t' -v p="$P" 'NR==1{for(i=1;i<=NF;i++) if($i==p) c=i; next} $c==1{print $1}' \
    "$WDIR/sp_parts.txt" > "$WDIR/sp_keep_$P.txt"
  bin/adjutant --pfile "$PFILE_PFX" vzs \
    --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
    --extract "$WDIR/sp_keep_$P.txt" \
    --out "$WDIR/sp_ex_$P" > "$WDIR/sp_ex_${P}_cli.txt" 2>&1 || true
  need "$WDIR/sp_ex_$P.sscore" "--score-partitions: the --extract run for $P produced nothing"
  for M in W1 W2 W3; do
    spsame "$WDIR/sp_ex_$P.sscore" "${M}_AVG" "$WDIR/sp.sscore" "${M}_${P}_AVG" \
      || SP_PARTBAD="$SP_PARTBAD ${M}/${P}"
  done
done
if [ -z "$SP_PARTBAD" ]; then
  echo "PASS  --score-partitions: every partition equals --extract on its own variants"
else
  echo "FAIL  --score-partitions: partition differs from its --extract equivalent for:$SP_PARTBAD"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 83-86 corroborate rather than add independent coverage.
#
# Additivity, the all-ones identity and the composition all FOLLOW from cases
# 81 and 82: once every partition equals its --extract equivalent and the total
# equals plain --score, these cannot be false.  Mutation confirms it -- every
# perturbation that breaks one of them trips 81 or 82 first, and the only
# mutation that leaves 81 and 82 intact is a uniform rescaling, which cancels in
# every relative comparison here and is score_check.sh's job to catch.
#
# They are kept because they state the properties a reader actually cares about
# in the form a reader would check them, and because they are nearly free.  But
# they should not be counted as independent evidence.
#
# Case 83: partitions that tile the variant set have additive sums.
#
# PHALF and POTHER are complements, so their _SUMs must add to the total _SUM.
# Sums, not averages: with per-partition denominators the averages deliberately
# do not add, and asserting on _AVG here would be asserting the wrong thing.
# ---------------------------------------------------------------------------
SP_ADD=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="W1_SUM")t=i; if($i=="W1_PHALF_SUM")a=i; if($i=="W1_POTHER_SUM")b=i}
                            if(!t||!a||!b){print "missing columns"; exit} next}
  { d=($a+$b)-$t; if(d<0) d=-d; if(d>1e-9) bad++ }
  END{ print bad+0 }' "$WDIR/sp.sscore")
# And the control: a non-tiling pair must NOT add up, or the tolerance is doing
# the work rather than the arithmetic.
SP_NOADD=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="W1_SUM")t=i; if($i=="W1_PHALF_SUM")a=i; if($i=="W1_PSMALL_SUM")c=i}; next}
  { d=($a+$c)-$t; if(d<0) d=-d; if(d>1e-9) bad++ }
  END{ print bad+0 }' "$WDIR/sp.sscore")
if [ "$SP_ADD" = "0" ] && [ "$SP_NOADD" -gt 100 ]; then
  echo "PASS  --score-partitions: tiling partitions have additive sums ($SP_NOADD rows differ for a non-tiling pair)"
else
  echo "FAIL  --score-partitions: additivity wrong (tiling mismatches=$SP_ADD, non-tiling mismatches=$SP_NOADD)"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 84: an all-ones partition is the total.
#
# The total is implemented as a synthetic all-ones partition, so this checks
# that the synthetic one is not special-cased into something different from
# what a user-supplied equivalent produces.
# ---------------------------------------------------------------------------
awk 'NR==1{print "SNP\tALLONES"; next} {print $1"\t1"}' "$WDIR/w3.txt" > "$WDIR/sp_all.txt"
bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --score-partitions "$WDIR/sp_all.txt" \
  --out "$WDIR/sp_ones" > "$WDIR/sp_ones_cli.txt" 2>&1 || true
need "$WDIR/sp_ones.sscore" "--score-partitions: the all-ones run produced nothing"
if spsame "$WDIR/sp_ones.sscore" W1_AVG "$WDIR/sp_ones.sscore" W1_ALLONES_AVG &&
   spsame "$WDIR/sp_ones.sscore" W1_AVG "$WDIR/sp_plain.sscore" W1_AVG; then
  echo "PASS  --score-partitions: an all-ones partition equals the total, which equals plain --score"
else
  echo "FAIL  --score-partitions: all-ones partition is not the total"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 85: a variant absent from the mask contributes to the total only.
#
# The reason the total column earns its place: without it, variants outside
# every partition would leave no trace in the output at all.
# ---------------------------------------------------------------------------
awk 'NR==1{print "SNP\tPONE"; next} NR<=200{print $1"\t1"}' "$WDIR/w3.txt" > "$WDIR/sp_partial.txt"
bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3-5 \
  --score-partitions "$WDIR/sp_partial.txt" \
  --out "$WDIR/sp_partial" > "$WDIR/sp_partial_cli.txt" 2>&1 || true
need "$WDIR/sp_partial.sscore" "--score-partitions: the partial-cover run produced nothing"
SP_SAME=$(spsame "$WDIR/sp_partial.sscore" W1_AVG "$WDIR/sp_partial.sscore" W1_PONE_AVG && echo same || echo differ)
if [ "$SP_SAME" = "differ" ] &&
   spsame "$WDIR/sp_partial.sscore" W1_AVG "$WDIR/sp_plain.sscore" W1_AVG; then
  echo "PASS  --score-partitions: uncovered variants reach the total and not the partition"
else
  echo "FAIL  --score-partitions: partial cover -- total should equal plain --score and differ from the partition"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 86: the evaluation flags see partitions as ordinary models.
#
# This is what the single-matrix output buys, and it is the reason the feature
# does not write one file per partition the way --q-score-range does.
# ---------------------------------------------------------------------------
awk -F'\t' 'NR==1{print "#IID\tSPTRAIT"} NR>1{print $1"\t"(NR%97)/13.0}' "$WDIR/sp.sscore" > "$WDIR/sp_pheno.txt"
bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read --score-col-nums 3 \
  --score-partitions "$WDIR/sp_parts.txt" \
  --pheno "$WDIR/sp_pheno.txt" --pheno-name SPTRAIT \
  --give-best-pgs --give-corr-matrix \
  --out "$WDIR/sp_comp" > "$WDIR/sp_comp_cli.txt" 2>&1 || true
need "$WDIR/sp_comp.best_pgs" "--score-partitions: the composed run produced no .best_pgs"
need "$WDIR/sp_comp.pgs_corr" "--score-partitions: the composed run produced no .pgs_corr"
SP_RANKED=$(awk -F'\t' -v c="$C_PGS_ID" 'NR>1{print $c}' "$WDIR/sp_comp.best_pgs" | sort | tr '\n' ' ')
SP_CORRN=$(awk -F'\t' 'NR==1{print NF-1; exit}' "$WDIR/sp_comp.pgs_corr")
if [ "$SP_RANKED" = "W1 W1_PHALF W1_POTHER W1_PSMALL " ] && [ "$SP_CORRN" = "4" ]; then
  echo "PASS  --score-partitions: --give-best-pgs ranks partitions and --give-corr-matrix spans them"
else
  echo "FAIL  --score-partitions: composition wrong (ranked '$SP_RANKED', corr width $SP_CORRN)"; exit 1
fi

# ---------------------------------------------------------------------------
# Cases 87-89: the guards.
# ---------------------------------------------------------------------------
sp_refuse() {  # sp_refuse <label> <message fragment> <args...>
  _label=$1; _want=$2; shift 2
  if bin/adjutant "$@" > "$WDIR/sp_rej.txt" 2>&1; then
    echo "FAIL  --score-partitions: $_label -- accepted, expected refusal"; exit 1
  elif saidflat "$WDIR/sp_rej.txt" "$_want"; then
    echo "PASS  --score-partitions: refused $_label"
  else
    echo "FAIL  --score-partitions: $_label -- refused, but not with '$_want':"
    cat "$WDIR/sp_rej.txt"; exit 1
  fi
}
printf 'SNP\n' > "$WDIR/sp_nocols.txt"
awk 'NR==1{print "SNP\tPA\tPB"; next} {print $1"\t"(NR%2)"\t2"}' "$WDIR/w3.txt" > "$WDIR/sp_badval.txt"
awk 'NR==1{print "SNP\tPA\tPA"; next} {print $1"\t1\t0"}' "$WDIR/w3.txt" > "$WDIR/sp_dup.txt"
awk 'NR==1{print "SNP\tPEMPTY\tPB"; next} {print $1"\t0\t1"}' "$WDIR/w3.txt" > "$WDIR/sp_empty.txt"
awk 'NR==1{print "SNP\tPA\tPB"; next} NR==2{print $1"\t1"; next} {print $1"\t1\t0"}' "$WDIR/w3.txt" > "$WDIR/sp_short.txt"

sp_refuse "a file with no partition columns" "no partition columns" \
  --pfile "$PFILE_PFX" vzs --score "$WDIR/w3.txt" 1 2 header-read \
  --score-partitions "$WDIR/sp_nocols.txt" --out "$WDIR/sp_g1"
sp_refuse "a value that is not 0 or 1, naming it" "is not 0 or 1" \
  --pfile "$PFILE_PFX" vzs --score "$WDIR/w3.txt" 1 2 header-read \
  --score-partitions "$WDIR/sp_badval.txt" --out "$WDIR/sp_g2"
sp_refuse "a duplicate partition name, naming it" "duplicate partition name 'PA'" \
  --pfile "$PFILE_PFX" vzs --score "$WDIR/w3.txt" 1 2 header-read \
  --score-partitions "$WDIR/sp_dup.txt" --out "$WDIR/sp_g3"
sp_refuse "a partition that selects nothing, naming it" "'PEMPTY' contains no variant" \
  --pfile "$PFILE_PFX" vzs --score "$WDIR/w3.txt" 1 2 header-read \
  --score-partitions "$WDIR/sp_empty.txt" --out "$WDIR/sp_g4"
sp_refuse "a short data row, naming the line" "line 2 of .* has 1 partition column; the header names 2" \
  --pfile "$PFILE_PFX" vzs --score "$WDIR/w3.txt" 1 2 header-read \
  --score-partitions "$WDIR/sp_short.txt" --out "$WDIR/sp_g5"
sp_refuse "use without --score" "must be used with --score" \
  --pfile "$PFILE_PFX" vzs --score-partitions "$WDIR/sp_parts.txt" --out "$WDIR/sp_g6"

# ---------------------------------------------------------------------------
# Cases 90-95: PGS IDs stay distinct when _AVG and _SUM name the same model.
#
# A PGS ID is the score column name minus plink2's _AVG / _SUM suffix.  With
# cols=+scoresums -- which writes BOTH columns for every model -- stripping
# unconditionally collapsed two columns onto one ID.  .pgs_corr came out
# headed `#ID W1 W2 W1 W2`, with the same names on its rows: ambiguous to any
# reader indexing by name, and R's read.table(row.names=1) refuses it outright.
# .best_pgs carried two rows per (phenotype, model) that looked identical and
# held different numbers, and --pgs-pheno-pairs matched both.
#
# The rule now: a column keeps its suffix when dropping it would make its ID
# equal to any other column's ID *or* to any other column's full name.  Both
# routes go through one implementation, per the two-routes-one-validator rule.
#
# The other half of the rule -- strip when the shortened name IS unique -- is
# what Case 0b ($id=="W4") and Case 70 (W1/PQT1 ...) already assert, so it is
# not restated here.  Break it and those go red, not these.
# ---------------------------------------------------------------------------
bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read cols=+scoresums --score-col-nums 3-5 \
  --out "$WDIR/dup" > "$WDIR/dup_cli.txt" 2>&1 || true
need "$WDIR/dup.sscore" "PGS ID collision: the cols=+scoresums run produced nothing"

# The .pgs_corr column names (header minus the leading #ID) and row names, each
# flattened to one line.  Kept separate because the bug put duplicates in both
# and a check on one alone would have passed while the other stayed broken.
corrcols() { head -1 "$1" | cut -f2- | tr '\t' '\n' | tr '\n' ' '; }
corrrows() { awk -F'\t' 'NR>1{print $1}' "$1" | tr '\n' ' '; }

# ---------------------------------------------------------------------------
# Case 90: --sscore route -- both flavours survive as their own model.
# ---------------------------------------------------------------------------
bin/adjutant --sscore "$WDIR/dup.sscore" \
  --give-corr-matrix --out "$WDIR/dup_ss" > "$WDIR/dup_ss_cli.txt" 2>&1 || true
need "$WDIR/dup_ss.pgs_corr" "PGS ID collision: the --sscore run produced no matrix"

DUP_WANT='W1_AVG W2_AVG W3_AVG W1_SUM W2_SUM W3_SUM '
DUP_SS_COLS=$(corrcols "$WDIR/dup_ss.pgs_corr")
DUP_SS_ROWS=$(corrrows "$WDIR/dup_ss.pgs_corr")
if [ "$DUP_SS_COLS" = "$DUP_WANT" ] && [ "$DUP_SS_ROWS" = "$DUP_WANT" ]; then
  echo "PASS  PGS ID collision: --sscore .pgs_corr names every column and row distinctly"
else
  echo "FAIL  PGS ID collision: --sscore .pgs_corr names wrong"
  echo "  expected: $DUP_WANT"
  echo "  cols:     $DUP_SS_COLS"
  echo "  rows:     $DUP_SS_ROWS"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 91: --pfile route, same file, same names.
#
# Not a duplicate of Case 90: --pfile reads its own .sscore back through a
# different function, and that copy had the same collapse.  Case 67 pins the
# two routes to the same matrix; this pins them to the same names.
# ---------------------------------------------------------------------------
bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w3.txt" 1 2 header-read cols=+scoresums --score-col-nums 3-5 \
  --give-corr-matrix --out "$WDIR/dup_pf" > "$WDIR/dup_pf_cli.txt" 2>&1 || true
need "$WDIR/dup_pf.pgs_corr" "PGS ID collision: the --pfile run produced no matrix"

DUP_PF_COLS=$(corrcols "$WDIR/dup_pf.pgs_corr")
DUP_PF_ROWS=$(corrrows "$WDIR/dup_pf.pgs_corr")
if [ "$DUP_PF_COLS" = "$DUP_WANT" ] && [ "$DUP_PF_ROWS" = "$DUP_WANT" ]; then
  echo "PASS  PGS ID collision: --pfile .pgs_corr names every column and row distinctly"
else
  echo "FAIL  PGS ID collision: --pfile .pgs_corr names wrong"
  echo "  expected: $DUP_WANT"
  echo "  cols:     $DUP_PF_COLS"
  echo "  rows:     $DUP_PF_ROWS"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 92: .best_pgs rows are addressable -- one per (PGS_ID, PHENO).
# ---------------------------------------------------------------------------
bin/adjutant --sscore "$WDIR/dup.sscore" \
  --pheno "$WDIR/pheno.txt" --pheno-name QT1 \
  --give-best-pgs --out "$WDIR/dup_bp" > "$WDIR/dup_bp_cli.txt" 2>&1 || true
need "$WDIR/dup_bp.best_pgs" "PGS ID collision: the .best_pgs run produced nothing"

DUP_BP_SET=$(pairset "$WDIR/dup_bp.best_pgs")
DUP_BP_EXPECT='W1_AVG/QT1 W1_SUM/QT1 W2_AVG/QT1 W2_SUM/QT1 W3_AVG/QT1 W3_SUM/QT1 '
if [ "$DUP_BP_SET" = "$DUP_BP_EXPECT" ]; then
  echo "PASS  PGS ID collision: .best_pgs has one row per (PGS_ID, PHENO)"
else
  echo "FAIL  PGS ID collision: .best_pgs pair set wrong"
  echo "  expected: $DUP_BP_EXPECT"
  echo "  got:      $DUP_BP_SET"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 93: --pgs-pheno-pairs selects one flavour, not both.
#
# The consequence the duplicate IDs had downstream: naming W1 matched the _AVG
# and _SUM rows alike, so a pair list could not ask for one of them.
#
# No --pheno-name: the pairs file already selects the phenotype and passing
# both is refused.  pheno.txt holds QT1 alone, so the pair list is the only
# thing narrowing anything here.
# ---------------------------------------------------------------------------
printf 'QT1\tW2_SUM\n' > "$WDIR/dup_pairs.txt"
bin/adjutant --sscore "$WDIR/dup.sscore" \
  --pheno "$WDIR/pheno.txt" \
  --give-best-pgs --pgs-pheno-pairs "$WDIR/dup_pairs.txt" \
  --out "$WDIR/dup_pp" > "$WDIR/dup_pp_cli.txt" 2>&1 || true
need "$WDIR/dup_pp.best_pgs" "PGS ID collision: the pairs-filtered run produced nothing"

DUP_PP_SET=$(pairset "$WDIR/dup_pp.best_pgs")
if [ "$DUP_PP_SET" = 'W2_SUM/QT1 ' ]; then
  echo "PASS  PGS ID collision: --pgs-pheno-pairs picks one flavour of a model"
else
  echo "FAIL  PGS ID collision: --pgs-pheno-pairs matched '$DUP_PP_SET', wanted 'W2_SUM/QT1 '"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 94: only the columns that actually collide keep their suffix.
#
# --score-col-nums can select a subset in which some models have both flavours
# and some have one.  W2 is alone here, so it keeps the short ID; W1 is not.
# This is what makes the rule per-column rather than all-or-nothing, and it is
# the only case that would notice the whole-set rule instead.
#
# Column numbers are literal indices into dup.sscore, whose fixed block is
# #IID SuperPop Population ALLELE_CT NAMED_ALLELE_DOSAGE_SUM: 6=W1_AVG,
# 7=W2_AVG, 9=W1_SUM.
# ---------------------------------------------------------------------------
bin/adjutant --sscore "$WDIR/dup.sscore" --score-col-nums 6,7,9 \
  --give-corr-matrix --out "$WDIR/dup_mix" > "$WDIR/dup_mix_cli.txt" 2>&1 || true
need "$WDIR/dup_mix.pgs_corr" "PGS ID collision: the mixed-subset run produced no matrix"

DUP_MIX_COLS=$(corrcols "$WDIR/dup_mix.pgs_corr")
if [ "$DUP_MIX_COLS" = 'W1_AVG W2 W1_SUM ' ]; then
  echo "PASS  PGS ID collision: a model with one flavour still gets the short ID"
else
  echo "FAIL  PGS ID collision: mixed subset named '$DUP_MIX_COLS', wanted 'W1_AVG W2 W1_SUM '"; exit 1
fi

# ---------------------------------------------------------------------------
# Case 95: a stripped name must not land on another column's FULL name.
#
# Weights columns W1 and W1_AVG make a .sscore holding W1_AVG, W1_AVG_AVG,
# W1_SUM, W1_AVG_SUM (6-9).  Selecting 6-8 leaves W1_AVG_AVG the only member of
# its pair, so the collide-with-another-ID clause alone would strip it to
# W1_AVG -- straight onto the W1_AVG that just kept its own suffix.  This is
# the case that makes the second clause do something.
# ---------------------------------------------------------------------------
awk 'BEGIN{print "ID\tA1\tW1\tW1_AVG"}
     !/^#/{print $3"\t"$5"\t"(NR%7)/10"\t"(NR%5)/10}' \
  Tests/TEST_PHASED_VCF/mt_mixed_phase.pvar | head -400 > "$WDIR/w_odd.txt"

bin/adjutant --pfile "$PFILE_PFX" vzs \
  --score "$WDIR/w_odd.txt" 1 2 header-read cols=+scoresums --score-col-nums 3-4 \
  --out "$WDIR/odd" > "$WDIR/odd_cli.txt" 2>&1 || true
need "$WDIR/odd.sscore" "PGS ID collision: the W1/W1_AVG run produced nothing"

ODD_COLS=$(head -1 "$WDIR/odd.sscore" | cut -f6-9 | tr '\t' '\n' | tr '\n' ' ')
if [ "$ODD_COLS" = 'W1_AVG W1_AVG_AVG W1_SUM W1_AVG_SUM ' ]; then
  echo "PASS  PGS ID collision: the W1/W1_AVG fixture has the columns this case needs"
else
  echo "FAIL  PGS ID collision: W1/W1_AVG fixture columns are '$ODD_COLS'"; exit 1
fi

bin/adjutant --sscore "$WDIR/odd.sscore" --score-col-nums 6-8 \
  --give-corr-matrix --out "$WDIR/odd_ss" > "$WDIR/odd_ss_cli.txt" 2>&1 || true
need "$WDIR/odd_ss.pgs_corr" "PGS ID collision: the W1/W1_AVG subset produced no matrix"

ODD_IDS=$(corrcols "$WDIR/odd_ss.pgs_corr")
if [ "$ODD_IDS" = 'W1_AVG W1_AVG_AVG W1_SUM ' ]; then
  echo "PASS  PGS ID collision: a stripped name never lands on another column's full name"
else
  echo "FAIL  PGS ID collision: named '$ODD_IDS', wanted 'W1_AVG W1_AVG_AVG W1_SUM '"; exit 1
fi
