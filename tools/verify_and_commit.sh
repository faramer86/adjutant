#!/bin/sh
# Clean-build, run all four suites, commit the working tree on green and
# throw it away on red.
#
#   tools/verify_and_commit.sh "prune(foo.cc): drop the bar cluster"
#
# The clean build is not optional.  An incremental build leaves stale .o files
# for every TU that includes a pruned header, so the link resolves symbols that
# no longer exist and the suites pass against a binary that cannot be rebuilt.
set -e
MSG=$1
[ -n "$MSG" ] || { echo "usage: $0 <commit message>"; exit 1; }

git diff --quiet && git diff --cached --quiet && { echo "SKIP  nothing to verify"; exit 0; }

# On red this throws the working tree away, so an uncommitted fix to the tool
# that produced the kill list would go with it -- which is exactly the edit you
# want to keep when a run fails. Commit tools/ first.
git diff --quiet -- tools/ || {
  echo "REFUSING  tools/ has uncommitted changes; commit them before verifying"
  exit 1
}

make clean > /dev/null 2>&1
if make -j8 BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD" \
        LDFLAGS="-ldl -lpthread -lz -flto" adjutant > /tmp/vac_build.log 2>&1 \
   && SKIP_BUILD=1 ./score_check.sh > /tmp/vac_score.log 2>&1 \
   && SKIP_BUILD=1 ./pgs_eval_check.sh > /tmp/vac_pgs_eval.log 2>&1 \
   && ./input_matrix.sh > /tmp/vac_matrix.log 2>&1 \
   && ./demo_give_best_pgs.sh > /tmp/vac_demo.log 2>&1; then
  git add -A
  git commit -q -m "$MSG"
  echo "PASS  $MSG"
else
  echo "FAIL  reverting working tree"
  grep -E " error:|Undefined symbols|^  \"" /tmp/vac_build.log | head -6 | sed 's/^/        /'
  # Leading-space count differs per suite -- score_check writes "FAIL", the
  # matrix "  FAIL", the demo "   FAIL" -- so match any indent rather than
  # enumerate them, or the failing suite reports nothing and the run looks
  # like a build break.
  grep -E "^ *FAIL|failed [1-9]" /tmp/vac_score.log /tmp/vac_pgs_eval.log \
       /tmp/vac_matrix.log /tmp/vac_demo.log 2>/dev/null | head -4 | sed 's/^/        /'
  git checkout -- .
  exit 1
fi
