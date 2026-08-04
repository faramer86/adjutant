#!/bin/sh
# Strip unreferenced typedefs/flagsets from one file, clean-build, run both
# suites, commit on green, revert on red.
set -e
F=$1
BEFORE=$(wc -l < "$F" | tr -d ' ')
python3 tools/deadtype.py strip "$F"
AFTER=$(wc -l < "$F" | tr -d ' ')
[ "$BEFORE" = "$AFTER" ] && { echo "SKIP  $F"; exit 0; }
make clean > /dev/null 2>&1
if make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
        LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /tmp/pt.log 2>&1 \
   && SKIP_BUILD=1 ./score_check.sh > /tmp/pts.log 2>&1 \
   && ./input_matrix.sh > /tmp/ptm.log 2>&1; then
  git add -A
  git commit -q -m "prune($F): drop unreferenced types, $BEFORE -> $AFTER lines"
  echo "PASS  $F  $BEFORE -> $AFTER lines  (committed)"
else
  echo "FAIL  $F  reverting"
  grep -E " error:" /tmp/pt.log | head -2 | sed 's/^/        /'
  git checkout -- "$F"; exit 1
fi
