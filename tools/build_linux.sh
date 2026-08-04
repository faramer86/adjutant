#!/bin/sh
# Build and verify an adjutant Linux binary inside a Debian container.
#
#   build_linux.sh <docker-platform> <output-filename>
#
# Debian bullseye is the base on purpose: glibc 2.31 (2020) is old enough that
# the result runs on the long-lived cluster images this tool is aimed at, while
# a newer base would raise that floor for no benefit.
#
# The binary is verified in the container that built it -- score_check.sh needs
# only in-repo Tests/ data, so shipping a compile-only artefact that was never
# run is avoidable and therefore not acceptable.
set -e
PLATFORM=$1
NAME=$2                       # output filename, written into <repo>/dist/

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
OUT=$SRC/dist                 # must live under $HOME: colima shares only those
mkdir -p "$OUT"               # paths, and a /tmp mount arrives silently empty

# The comparison helper is passed through the environment rather than a bind
# mount: colima only shares paths under $HOME, so a /tmp mount silently arrives
# empty in the container.
HELPER_B64=$(base64 < "$SRC/tools/sscore_cmp.sh" | tr -d '\n')

docker run --rm --platform "$PLATFORM" \
  -v "$SRC":/src:ro -v "$OUT":/out \
  -e HELPER_B64="$HELPER_B64" -e OUTNAME="$NAME" \
  debian:bullseye bash -euc '
    export DEBIAN_FRONTEND=noninteractive
    apt-get update -qq
    apt-get install -y -qq --no-install-recommends \
      g++ make zlib1g-dev libopenblas-dev liblapacke-dev file > /dev/null

    # Copy out of the read-only mount: the build writes .o files beside the
    # sources, and the host tree holds macOS objects we must not disturb.
    mkdir -p /build && cp -r /src/. /build/ && cd /build
    rm -rf .git
    make clean > /dev/null 2>&1 || true

    # -static-libstdc++/-static-libgcc so the binary does not carry a dependency
    # on the build image toolchain runtime.  BLAS/LAPACK and glibc stay dynamic.
    make -j"$(nproc)" \
      BASEFLAGS="-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD" \
      BLASFLAGS64="-lopenblas -llapacke" \
      LDFLAGS="-ldl -lpthread -lz -lm -static-libstdc++ -static-libgcc" \
      adjutant > /tmp/build.log 2>&1

    echo "=== built ==="
    file bin/adjutant
    ./bin/adjutant --version
    echo "=== runtime dependencies ==="
    ldd bin/adjutant | sed "s/^/  /"

    # score_check.sh diffs .log/.stdout as well, and score_regression.sh
    # normalises only the volatile fields macOS emits.  Linux adds
    # ", ~NNNN available" to the RAM line and a "Using up to N compute threads."
    # line, so the suite fails on log formatting even when every number matches.
    # Run it for the record; let the .sscore comparison carry the verdict.
    echo "=== score_check.sh (log-format noise expected on Linux) ==="
    SKIP_BUILD=1 ./score_check.sh ./bin/adjutant > /tmp/sc.log 2>&1 || true
    grep -E "^(PASS|FAIL)" /tmp/sc.log | head -5 | sed "s/^/  /"

    echo "=== .sscore comparison (the numerical contract) ==="
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
