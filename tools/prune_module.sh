#!/bin/sh
# Prune one translation unit down to what is reachable from the symbols other
# objects actually reference, then prove the tool still works.  Reverts on any
# failure, commits on success.
#
#   tools/prune_module.sh <module-basename>       e.g. adjutant_common
set -e
M=$1
[ -f "$M.cc" ] || { echo "no such module: $M.cc"; exit 1; }

BEFORE=$(wc -l < "$M.cc" | tr -d ' ')
make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
     LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /dev/null 2>&1 || true
ROOTS=$(python3 tools/extroots.py "$M" $(ls *.o include/*.o 2>/dev/null))

HDRARG=""
[ -f "$M.h" ] && HDRARG="--header $M.h"
# shellcheck disable=SC2086
python3 tools/deadcode.py strip "$M.cc" "$ROOTS" $HDRARG

make clean > /dev/null 2>&1
if make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
        LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /tmp/prune_build.log 2>&1 \
   && SKIP_BUILD=1 ./score_check.sh > /tmp/prune_score.log 2>&1 \
   && ./input_matrix.sh > /tmp/prune_matrix.log 2>&1; then
  AFTER=$(wc -l < "$M.cc" | tr -d ' ')
  git add -A
  git commit -q -m "prune($M): drop unreachable functions, $BEFORE -> $AFTER lines"
  echo "PASS  $M.cc  $BEFORE -> $AFTER lines  (committed)"
else
  echo "FAIL  $M.cc  reverting"
  grep -E " error:|Undefined symbols" /tmp/prune_build.log | head -3 | sed 's/^/        /'
  grep -E "^FAIL" /tmp/prune_score.log /tmp/prune_matrix.log 2>/dev/null | head -3 | sed 's/^/        /'
  git checkout -- "$M.cc"
  [ -f "$M.h" ] && git checkout -- "$M.h"
  exit 1
fi
