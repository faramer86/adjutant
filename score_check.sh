#!/bin/sh
# Verify an adjutant build still behaves exactly like pristine plink2 for
# --score, and that the removed commands are genuinely gone.
#
# The contract this enforces, which must hold for every build:
#   1. --pfile  (.pgen/.pvar/.psam)  loads and scores
#   2. --bfile  (.bed/.bim/.fam)     loads and scores
#   3. --score  with every modifier, plus --score-col-nums
#   and all three produce output byte-identical to stock plink2.
#
# Usage:  ./score_check.sh [adjutant-binary]
#
# Exit 0 = pass.

set -e

BIN=${1:-./bin/adjutant}
GOLDEN=score_regression_golden
OUT=${OUT_DIR:-/tmp/score_after}

# Rebuild first. Without this, a failed build silently leaves the previous
# binary in place and the whole check passes against stale output.
if [ "${SKIP_BUILD:-0}" != "1" ]; then
  if ! make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
       LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /tmp/score_build.log 2>&1; then
    echo "FAIL  build failed:"
    grep -E "error:|Undefined" /tmp/score_build.log | head -20
    exit 1
  fi
fi

./score_regression.sh "$BIN" "$OUT" > /dev/null 2>&1

fail=0

# 1. Every kept case must be byte-identical to the pristine baseline.
#    gone_* cases are excluded: they are *supposed* to change.
if diff -r -x 'gone_*' -x '_failures.txt' "$GOLDEN" "$OUT" > /tmp/score_check.diff 2>&1; then
  kept=$(find "$OUT" -type f ! -name 'gone_*' ! -name '_failures.txt' | wc -l | tr -d ' ')
  echo "PASS  --score behaviour identical to baseline ($kept files)"
else
  echo "FAIL  --score behaviour diverged from baseline:"
  head -40 /tmp/score_check.diff
  fail=1
fi

# 2. Every removed command must be rejected, not silently accepted.
for f in "$OUT"/gone_*.stdout; do
  [ -e "$f" ] || continue
  name=$(basename "$f" .stdout)
  if grep -qE 'Error: Unrecognized flag|Error: Invalid flag' "$f"; then
    :
  else
    echo "FAIL  $name was not rejected as an unrecognized flag:"
    tail -3 "$f" | sed 's/^/        /'
    fail=1
  fi
done
[ "$fail" -eq 0 ] && echo "PASS  all removed commands rejected as unrecognized flags"

# Spell out the contract, so a green run says what it actually proved.
if [ "$fail" -eq 0 ]; then
  pf=$(ls "$OUT" | grep -c '^basic\.sscore$' || true)
  bf=$(ls "$OUT" | grep -c '^bed_basic\.sscore$' || true)
  if [ "$pf" -eq 1 ] && [ "$bf" -eq 1 ]; then
    echo "PASS  contract: --pfile, --bfile, --score + all modifiers"
  else
    echo "FAIL  contract: an input format produced no .sscore (pfile=$pf bfile=$bf)"
    exit 1
  fi
fi

exit $fail
