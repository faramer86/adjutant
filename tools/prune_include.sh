#!/bin/sh
# Same workflow as prune_module.sh, for a translation unit under include/.
set -e
M=$1
[ -f "include/$M.cc" ] || { echo "no such module: include/$M.cc"; exit 1; }
BEFORE=$(wc -l < "include/$M.cc" | tr -d ' ')
make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
     LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /dev/null 2>&1 || true
ROOTS=$(python3 tools/extroots.py "include/$M" $(ls *.o include/*.o 2>/dev/null))
HDRARG=""
[ -f "include/$M.h" ] && HDRARG="--header include/$M.h"
# shellcheck disable=SC2086
python3 tools/deadcode.py strip "include/$M.cc" "$ROOTS" $HDRARG
make clean > /dev/null 2>&1
if make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
        LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /tmp/pi_build.log 2>&1 \
   && SKIP_BUILD=1 ./score_check.sh > /tmp/pi_score.log 2>&1 \
   && ./input_matrix.sh > /tmp/pi_matrix.log 2>&1; then
  AFTER=$(wc -l < "include/$M.cc" | tr -d ' ')
  git add -A
  git commit -q -m "prune(include/$M): drop unreachable functions, $BEFORE -> $AFTER lines"
  echo "PASS  include/$M.cc  $BEFORE -> $AFTER lines  (committed)"
else
  echo "FAIL  include/$M.cc  reverting"
  grep -E " error:|Undefined symbols" /tmp/pi_build.log | head -2 | sed 's/^/        /'
  git checkout -- "include/$M.cc"
  [ -f "include/$M.h" ] && git checkout -- "include/$M.h"
  exit 1
fi
