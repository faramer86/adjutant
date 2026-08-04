#!/bin/sh
# Compare every .sscore against its golden and, where they differ, report how
# big the differing numbers actually are.
#
# "Not byte-identical" and "numerically wrong" are different claims.  A score of
# 1e-18 where the true value is 0 is cancellation residue; a difference in a
# meaningful digit is a bug.  Printing the magnitude is what separates them.
GOLD=${1:-score_regression_golden}
OUT=${2:-/tmp/score_after}
same=0; diffr=0; worst=0
for f in "$GOLD"/*.sscore; do
  b=$(basename "$f")
  [ -f "$OUT/$b" ] || continue
  if cmp -s "$f" "$OUT/$b"; then same=$((same+1)); continue; fi
  diffr=$((diffr+1))
  paste "$f" "$OUT/$b" | awk -F'\t' -v n="$b" '
    NR==1 { nf=NF/2; next }
    {
      for (i=1;i<=nf;i++) {
        x=$i; y=$(i+nf)
        if (x!=y) {
          ax=(x<0?-x:x); ay=(y<0?-y:y); m=(ax>ay?ax:ay)
          if (m>mx) mx=m
          c++
        }
      }
    }
    END { printf "    %-26s %4d differing fields, largest magnitude %g\n", n, c, mx }'
done
echo "  identical .sscore: $same   differing: $diffr"
