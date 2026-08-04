#!/bin/sh
# Strip unreferenced inline definitions from one header, then prove the tool
# still builds clean and passes both suites.  Reverts on failure.
set -e
H=$1
[ -f "$H" ] || { echo "no such header: $H"; exit 1; }
BEFORE=$(wc -l < "$H" | tr -d ' ')
python3 tools/deadheader.py strip "$H"
AFTER=$(wc -l < "$H" | tr -d ' ')
if [ "$BEFORE" = "$AFTER" ]; then echo "SKIP  $H (nothing to remove)"; exit 0; fi
make clean > /dev/null 2>&1
if make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
        LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /tmp/ph_build.log 2>&1 \
   && SKIP_BUILD=1 ./score_check.sh > /tmp/ph_score.log 2>&1 \
   && ./input_matrix.sh > /tmp/ph_matrix.log 2>&1; then
  git add -A
  git commit -q -m "prune($H): drop unreferenced inline defs, $BEFORE -> $AFTER lines"
  echo "PASS  $H  $BEFORE -> $AFTER lines  (committed)"
else
  echo "FAIL  $H  reverting"
  grep -E " error:|Undefined symbols" /tmp/ph_build.log | head -2 | sed 's/^/        /'
  git checkout -- "$H"
  exit 1
fi
