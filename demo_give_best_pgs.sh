#!/bin/sh
# Demonstrate and self-validate --give-best-pgs.
#
# Simulates a phenotype driven by a KNOWN score, then checks the feature both
# ranks that score first and recovers the coefficient it was built with.  A
# demo that only shows the command running proves nothing; this one fails
# loudly if the maths is wrong.
#
#   ./demo_give_best_pgs.sh
set -e
cd "$(dirname "$0")"

P=Tests/TEST_PHASED_VCF/mt_mixed_phase
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT

TRUE_BETA=1.2       # coefficient we simulate with
DRIVER=MODEL_B      # the score that should win

# Resolve a column index by name, from the file's own header.  This demo read
# BETA from a literal $5 and R2_NULL from $8; when .best_pgs grew from 16
# columns to 21 both silently moved and the demo failed reporting BETA NA.  A
# name cannot point at the wrong column, and a name that is not in the header
# stops the script instead of yielding $"".  pgs_eval_check.sh carries the same
# six lines -- both scripts have to run standalone.
col() {  # col <file> <COLUMN NAME> -> 1-based index
  _idx=$(awk -F'\t' -v want="$2" \
    'NR==1{for(i=1;i<=NF;i++){n=$i; sub(/^#/,"",n); if(n==want){print i; exit}}}' \
    "$1")
  [ -n "$_idx" ] || { echo "   FAIL  $1 has no column named '$2'" >&2; exit 1; }
  printf '%s\n' "$_idx"
}

echo "1. Three candidate PGS models -----------------------------------------"
awk 'BEGIN{print "ID\tA1\tMODEL_A\tMODEL_B\tMODEL_C"}
     !/^#/{split($5,a,","); print $3"\t"a[1]"\t"((NR*7)%13-6)/10"\t"((NR*3)%11-5)/10"\t"((NR*5)%7-3)/10}' \
  $P.pvar | head -450 > "$D/weights.txt"
head -3 "$D/weights.txt" | sed 's/^/   /'

echo "2. Score them to get the N x M PGS matrix ------------------------------"
bin/adjutant --pfile $P vzs --score "$D/weights.txt" 1 2 header-read \
  --score-col-nums 3-5 --out "$D/pgs" > /dev/null 2>&1
head -1 "$D/pgs.sscore" | tr '\t' ' ' | sed 's/^/   /'

echo "3. Simulate a phenotype driven by $DRIVER -----------------------------"
# Standardise the driving score first: raw scores are ~1e-3, so an unscaled
# coefficient would be swamped by the noise term and the demo would silently
# rank the wrong model.
DRIVER_COL=$(col "$D/pgs.sscore" "${DRIVER}_AVG")
read MEAN SD <<EOF
$(awk -v c="$DRIVER_COL" 'NR>1{s+=$c; ss+=$c*$c; n++} END{m=s/n; printf "%.12f %.12f", m, sqrt(ss/n-m*m)}' "$D/pgs.sscore")
EOF
echo "   $DRIVER: mean=$MEAN sd=$SD"
echo "   TRAIT = $TRUE_BETA * z($DRIVER) + 0.05 * AGE + noise"
awk -v m="$MEAN" -v sd="$SD" -v b="$TRUE_BETA" -v c="$DRIVER_COL" \
    'BEGIN{OFS="\t"; srand(11); print "#IID","TRAIT"}
     NR>1{z=($c-m)/sd; age=(NR%53)+20;
          noise=(rand()+rand()+rand()+rand()-2);
          print $1, b*z + 0.05*age + noise}' "$D/pgs.sscore" > "$D/pheno.txt"

awk 'BEGIN{OFS="\t"; print "#IID","AGE","PC1"}
     NR>1{print $1,(NR%53)+20,((NR*17)%100)/100.0}' "$D/pgs.sscore" > "$D/covar.txt"

echo "4. Rank the models ----------------------------------------------------"
bin/adjutant --pfile $P vzs \
  --pheno "$D/pheno.txt" --pheno-name TRAIT \
  --covar "$D/covar.txt" \
  --score "$D/weights.txt" 1 2 header-read --score-col-nums 3-5 \
  --give-best-pgs --out "$D/result" > /dev/null 2>&1
column -t "$D/result.best_pgs" | sed 's/^/   /'

echo "5. Validate -----------------------------------------------------------"
BEST="$D/result.best_pgs"
TOP=$(awk -F'\t' -v c="$(col "$BEST" PGS_ID)" 'NR==2{print $c}' "$BEST")
GOT=$(awk -F'\t' -v c="$(col "$BEST" BETA)"   'NR==2{print $c}' "$BEST")

if [ "$TOP" != "$DRIVER" ]; then
  echo "   FAIL  top model is $TOP, expected $DRIVER"; exit 1
fi
echo "   ok    $DRIVER ranked first"

# BETA is per standard deviation, and the phenotype was built from the
# standardised score, so the estimate should land near TRUE_BETA.
if awk -v g="$GOT" -v t="$TRUE_BETA" 'BEGIN{d=g-t; if(d<0)d=-d; exit !(d<0.15)}'; then
  echo "   ok    BETA $GOT recovers the simulated $TRUE_BETA"
else
  echo "   FAIL  BETA $GOT is far from the simulated $TRUE_BETA"; exit 1
fi

# Every row shares one null model, so R2_NULL must be identical throughout.
if [ "$(awk -F'\t' -v c="$(col "$BEST" R2_NULL)" 'NR>1{print $c}' "$BEST" \
        | sort -u | wc -l | tr -d ' ')" != "1" ]; then
  echo "   FAIL  R2_NULL differs between rows -- sample sets diverged"; exit 1
fi
echo "   ok    R2_NULL identical on every row"

echo
echo "PASS  --give-best-pgs ranks correctly and recovers the true coefficient"
