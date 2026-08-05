# adjutant binaries

**This directory holds no binaries between releases** — only this file, which is
the record of how they are produced and what must be checked before any of them
is published. Stale binaries that look current are worse than none: they carry
an older output layout while the README documents the new one, and nothing about
the file says so. They are therefore built fresh at release time, all targets
from one commit, and deleted afterwards.

Everything is compiled with `-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD`.

> **Do not distribute a build older than `4a18d56`.** Everything before it has a
> stack buffer overflow in the `.best_pgs` writer, reachable from an input file:
> a `.sscore` whose score-column name runs past about 3.6 kB overflows a fixed
> 4096-byte buffer. The pre-fix binary segfaults on a 9000-character name; this
> one writes the row correctly.

`-ffp-contract=off` is not optional. It stops the compiler fusing multiply-add
in the scoring dot product, which is what keeps `-O3` bit-identical to `-O2`.

Two targets are supported: **macOS Apple Silicon** and **Linux x86-64**. Earlier
builds for macOS Intel, the universal binary and Linux aarch64 have been dropped.

| file | platform | approx. size |
|---|---|---|
| `adjutant-macos-arm64` | macOS Apple Silicon | 1.3 MB |
| `adjutant-linux-x86_64` | Linux x86-64, glibc ≥ 2.31 | 1.6 MB |
| `adjutant-linux-x86_64-static` | Linux x86-64, **no dependencies** | 28 MB |

**Build every target from the same commit.** Mixing them is how a release ships
binaries that disagree with each other: `.matched` gained `CHR`/`POS`/`A2`
columns after `dc4908a`, so a set assembled across that boundary would have one
binary writing a layout the others cannot read.

## What to verify before publishing

Every binary must be **run**, not merely compiled.

`score_check.sh` covers `--score` and nothing else, so everything built on top of
it was exercised separately: plain and partitioned scoring, `--give-best-pgs`
with and without `--folds`, `--give-corr-matrix`, `--pgs-pheno-pairs`,
`--pgs-model`, `--pgs-models-list`, all four MHC regions on both routes, and the
refusal of the `--extract`/`--exclude` interval modifiers. Sixteen outputs in
all, compared by checksum across platforms.

- **macOS arm64** — all four suites green: 89/89 goldens byte-identical, the
  input matrix, every `pgs_eval_check` assertion, and the demo.
- **Linux dynamic vs Linux static** — every output byte-identical. Measured at
  `dc4908a`: 16 of 16.
- **macOS vs Linux** — everything byte-identical except `--score` `.sscore`,
  where a handful of near-zero cells differ; see below. Measured at `dc4908a`:
  14 of 16 identical, `.best_pgs`, `.pgs_corr`, `.matched`, `.match_report` and
  every MHC output among them.

A build verified only by `score_check.sh` has never executed `--give-best-pgs`,
`--give-corr-matrix`, `--pgs-pheno-pairs`, `--score-partitions` or the MHC
filters even once. Repeat the full sweep after any rebuild.

## Byte-identity is per-toolchain

The contract is that `--score` output is byte-identical to stock plink2. That
holds exactly for the **native macOS arm64 clang build** — the toolchain the
goldens were generated with. It does not hold bit-for-bit on Linux, which uses
gcc and OpenBLAS rather than clang and Accelerate.

The disagreement is confined to cells whose value is **itself numerically zero**.
Measured across 300 samples: the largest difference anywhere is **8.33e-18**, and
the largest value involved in any differing cell is **5.55e-18**. Those are
columns where the score cancels and only rounding residue survives. Every
meaningful digit of every other field matches, and the difference is unchanged
from the previously shipped Linux build — it is a property of the toolchain, not
of any code change.

## Linux runtime requirements

`adjutant-linux-x86_64` is dynamically linked against OpenBLAS. One package is
enough — `libopenblas0` pulls in `libgfortran5`:

```sh
apt-get install libopenblas0        # Debian/Ubuntu
```

Verified by running it in a bare `debian:bullseye`: without the package it fails
with `libopenblas.so.0: cannot open shared object file`; with it alone it starts
and reports its version.

`adjutant-linux-x86_64-static` needs nothing at all — verified by running it in
`busybox`, which carries no glibc whatsoever. It costs size: 27.8 MB against
1.62 MB. That is OpenBLAS's per-CPU kernels, not debug info — `strip
--strip-unneeded` recovers only about 1 MB, so this build is left unstripped.

Built on Debian bullseye (glibc 2.31, 2020) so the binaries run on older
long-lived cluster images.

## Notable changes to check against

`.matched` is now `CHR POS ID A1 A2` followed by one weight column per model,
every identifying field taken from the matched `.bim` row. `--score` therefore
needs column 3 as the ID, 4 as the allele, and `--score-col-nums` for the
weights — even for a single model. A binary predating that change writes
`ID A1 <weights>` instead.

Since `f32fd02`:

- **`--pgs-model` / `--pgs-models-list` can join on variant ID** — a leading
  `id` keyword, for weight files whose positions belong to a different build.

- **`--exclude-mhc-hg19` / `--exclude-xmhc-hg19` / `--exclude-mhc-hg38` /
  `--exclude-xmhc-hg38`** — drop chromosome 6 MHC variants, on both the scoring
  and the `--bim`/`--pgs-model` route. Note that xMHC is **not** a superset of
  MHC: the two come from different published definitions and the GRC region ends
  79,933 bp past the xMHC end, so the extended flag can exclude *fewer* variants.
- **The `--extract`/`--exclude` interval modifiers are now refused.**
  `bed0`, `ibed0`, `bed1`, `ibed1` and `range` used to parse, appear in "Options
  in effect", and then filter nothing — an `--extract` meant to restrict a score
  to one region silently returned the genome-wide score. They are hard errors
  now. ID-based `--extract`/`--exclude` are unaffected.
- **`--pgs-models-list` no longer scales memory with model count.** It used to
  hold every model's parsed records until the end, at 227 bytes per source
  variant — about 11x the size of the input files. 249 models totalling 5.9 GB
  needed roughly 65 GB and were OOM-killed. Measured on 10M source variants:
  2.08 GB down to 0.30 GB, and peak now grows about 2 MB per model instead of 20.

Earlier additions still present: `--pgs-models-list`, `--give-corr-matrix`,
`--pgs-pheno-pairs`, `--score-partitions`, `--resolve-ambiguous` (renamed from
`--resolve-ambig`; the old spelling is rejected), and the `--sscore` fixes for
PLINK's `#FID IID` header and for identifying score columns by name.

## Reproducing

macOS, native:

```sh
make release                                     # -> bin/adjutant
```

Linux, in a Debian bullseye container:

```sh
# dynamic
make -j"$(nproc)" \
  BASEFLAGS="-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD" \
  BLASFLAGS64="-lopenblas -llapacke" \
  LDFLAGS="-ldl -lpthread -lz -lm -static-libstdc++ -static-libgcc" adjutant

# static -- link order matters, the linker resolves left to right
make -j"$(nproc)" \
  BASEFLAGS="-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD" \
  BLASFLAGS64="-llapacke -lopenblas -lgfortran -lquadmath" \
  LDFLAGS="-static -ldl -lpthread -lz -lm" adjutant
```

## Traps worth recording

**Never `cp` over a binary in place.** `cp` overwrites the existing inode, and
overwriting the text pages of an executable another process is still running
wedges that process in uninterruptible state on macOS. Every later exec of the
same inode then wedges too, at zero CPU, which looks exactly like an application
hang — and the wedged processes cannot be killed. `rm -f` first, so a running
process keeps the old inode and the new file gets a fresh one.

**`score_check.sh` cannot pass unmodified on Linux.** `score_regression.sh`
normalises only the volatile fields macOS emits; Linux adds `, ~NNNN available`
to the RAM line and a `Using up to N compute threads.` line. That is log
formatting, not scoring. The `.sscore` files are the numerical contract and are
compared directly.

**Cross-compiling to x86-64 on Apple Silicon is unreliable.** Under qemu,
`cc1`/`g++` take a SIGSEGV at random — a different translation unit each time,
on files that compile fine natively. `make` resumes from the object files
already produced, so retrying converges; drop to `-j1` after the first failure.
Verify the result from the ELF header (`7f454c46` magic, `e_machine` `0x3e`)
rather than with `file`, which is **not installed** in `debian:bullseye` — a
check that silently rejected 19 good builds. And exclude the host's `bin/`
from whatever you copy into the container: a stray macOS binary landing there
turns an existence check into a false pass.
