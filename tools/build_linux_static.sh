#!/bin/sh
# Fully static Linux build: no libopenblas, no libgfortran, no glibc at runtime.
#
#   build_linux_static.sh <docker-platform> <output-filename>
#
# Static archives for everything we need ship in debian:bullseye, so the only
# real work is link order -- with -static the linker resolves left to right, and
# LAPACKE must precede OpenBLAS, which must precede the Fortran runtime.
set -e
PLATFORM=$1
NAME=$2

# Both arguments are required.  Without them docker falls back to the HOST
# platform -- aarch64 on an Apple-silicon colima -- and builds an ARM binary
# that would land under an x86_64 filename.  Caught once, three minutes in, only
# because the copy failed on an empty OUTNAME; check up front instead.
if [ -z "$PLATFORM" ] || [ -z "$NAME" ]; then
  echo "usage: $0 <docker-platform> <output-filename>" >&2
  echo "  e.g. $0 linux/amd64 adjutant-linux-x86_64" >&2
  exit 2
fi
# Repo root, derived from this script's location: a hardcoded path here
# would work only on the machine it was written on.
SRC=$(cd "$(dirname "$0")/.." && pwd)
OUT=$SRC/dist                 # under $HOME: colima shares only those paths
mkdir -p "$OUT"

HELPER_B64=$(base64 < "$SRC/tools/sscore_cmp.sh" | tr -d '\n')

docker run --rm --platform "$PLATFORM" \
  -v "$SRC":/src:ro -v "$OUT":/out \
  -e HELPER_B64="$HELPER_B64" -e OUTNAME="$NAME" \
  debian:bullseye bash -euc '
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends \
      g++ make zlib1g-dev libopenblas-dev liblapacke-dev gfortran file > /dev/null

    mkdir -p /build && cp -r /src/. /build/ && cd /build
    rm -rf .git
    make clean > /dev/null 2>&1 || true

    # libquadmath exists on x86_64 but not aarch64; only ask for it where it is.
    QUAD=""
    ls /usr/lib/gcc/*/*/libquadmath.a > /dev/null 2>&1 && QUAD="-lquadmath"

    make -j"$(nproc)" \
      BASEFLAGS="-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD" \
      BLASFLAGS64="-llapacke -lopenblas -lgfortran $QUAD" \
      LDFLAGS="-static -ldl -lpthread -lz -lm" \
      adjutant > /tmp/build.log 2>&1 || { echo "BUILD FAILED"; tail -25 /tmp/build.log; exit 1; }

    echo "=== built ==="
    file bin/adjutant
    ./bin/adjutant --version
    echo "=== dynamic dependencies (want: none) ==="
    ldd bin/adjutant 2>&1 | sed "s/^/  /"

    echo "=== .sscore comparison (the numerical contract) ==="
    SKIP_BUILD=1 ./score_check.sh ./bin/adjutant > /tmp/sc.log 2>&1 || true
    echo "$HELPER_B64" | base64 -d > /tmp/cmp.sh
    sh /tmp/cmp.sh score_regression_golden /tmp/score_after

    # Unlink before copying.  cp overwrites IN PLACE, and overwriting the
    # text pages of an executable that another process is still running wedges
    # that process in uninterruptible state on macOS -- after which every
    # later exec of the same inode wedges too, at zero CPU, looking exactly
    # like an application hang.  Removing first leaves the running process on
    # the old inode and gives the new file a fresh one.
    rm -f "/out/$OUTNAME"
    cp bin/adjutant "/out/$OUTNAME"
    echo "=== done ==="
  '
