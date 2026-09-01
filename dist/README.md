# adjutant binaries

**Built binaries live here.** This file records how they are produced and what
must be checked before any is published.

**Rebuild every target from the same commit, together.** A stale binary that
looks current is worse than none — it writes an older output layout while the
README documents the new one, and nothing about the file says which. `.matched`
gained `CHR`/`POS`/`A2` columns after `dc4908a`, so a set assembled across that
boundary would have one binary emitting a layout the others cannot read.

Everything is compiled with `-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD`.

> **Do not distribute a build older than `4a18d56`.** Everything before it has a
> stack buffer overflow in the `.best_pgs` writer, reachable from an input file:
> a `.sscore` whose score-column name runs past about 3.6 kB overflows a fixed
> 4096-byte buffer. The pre-fix binary segfaults on a 9000-character name; this
> one writes the row correctly.

`-ffp-contract=off` is not optional. It stops the compiler fusing multiply-add
in the scoring dot product, which is what keeps `-O3` bit-identical to `-O2`.

Two targets ship, both needing nothing installed:

| file | platform | approx. size |
|---|---|---|
| `adjutant-macos-arm64` | macOS Apple Silicon | 1.3 MB |
| `adjutant-linux-x86_64-static` | Linux x86-64 | 28 MB |

Earlier builds for macOS Intel, the universal binary and Linux aarch64 have been
dropped, as has the dynamically linked Linux build — the static one runs
anywhere and needs no `libopenblas0`. Its recipe is still below if you want it.

### "Static" means different things on the two platforms

The Linux binary is genuinely `-static`: one file, no interpreter, verified in
`busybox`, which carries no glibc at all.

**A fully static macOS binary is not possible.** Apple ships no static
`libSystem` and has not supported statically linked executables for two decades,
so `-static` cannot link there. What the macOS build does instead is depend on
**nothing outside the OS**:

```
Accelerate.framework   libSystem.B.dylib   libz.1.dylib   libc++.1.dylib
```

All four are present on every macOS. There is nothing to install, which is the
property that matters for distribution — reached by a different route, and worth
stating precisely rather than calling the result "static".

## What to verify before publishing

Every binary must be **run**, not merely compiled.

**The commit gate does not build what ships.** `tools/verify_and_commit.sh`
compiles with `BASEFLAGS="-g -ffp-contract=off -DZSTD_MULTITHREAD"` — no `-O3`,
no `-DNDEBUG` — while `make release` adds both. A green gate is therefore
evidence about the unoptimised binary, and says nothing about the optimised one
beyond what they share. Run the suites and the sweep again against the actual
`make release` output before copying it here. The two are easy to tell apart by
size: the release build is 1,330,792 bytes, the gate's is 1,449,416.

`make release` is reproducible — two clean builds of `b2f71f8` gave the same
SHA-256 — so a hash mismatch against `SHA256SUMS` means a different commit or a
different recipe, not build nondeterminism.

`score_check.sh` covers `--score` and nothing else, so everything built on top of
it was exercised separately: plain and partitioned scoring, `--give-best-pgs`
with and without `no-na`, `--give-corr-matrix`, `--pgs-pheno-pairs`,
`--pgs-model` on both keys, `--pgs-models-list`, all four MHC regions on both
routes, a `.sscore` wider than the old 1024-column ceiling, and the refusal of
the `--extract`/`--exclude` interval modifiers. Sixteen commands, 27 outputs,
compared by checksum across platforms.

**Make the sweep assert on each step.** The first version of it ran each command
with `|| true` and checksummed whatever files appeared. Five of the sixteen were
failing on wrong flag names and column specs — `--fold-seed` does not exist, and
`--pgs-models-list` could not resolve the bare filenames in `models.txt` — and
the run still reported success, because a missing output is indistinguishable
from a skipped comparison. Name the file each step must produce and check it.

- **macOS arm64** — all four suites green: 89/89 goldens byte-identical, the
  input matrix, every `pgs_eval_check` assertion, and the demo.
- **macOS vs Linux** — everything byte-identical except `--score` `.sscore`,
  where a handful of near-zero cells differ; see below. `.matched`,
  `.match_report` (both key modes), the merged matrix and every MHC output all
  agree exactly.

A build verified only by `score_check.sh` has never executed `--give-best-pgs`,
`--give-corr-matrix`, `--pgs-pheno-pairs`, `--score-partitions` or the MHC
filters even once. Repeat the full sweep after any rebuild.

## Byte-identity is per-toolchain

The contract is that `--score` output is byte-identical to stock plink2. That
holds exactly for the **native macOS arm64 clang build** — the toolchain the
goldens were generated with. It does not hold bit-for-bit on Linux, which uses
gcc and OpenBLAS rather than clang and Accelerate.

The disagreement is confined to cells whose value is **itself numerically zero**.
Re-measured at every release, across the six `.sscore` outputs of the sweep:
**49 of 16,000 numeric cells** differ, by at most **4.16e-18**, and the largest
value appearing in any differing cell is **1.04e-17**. Those are columns where
the score cancels and only rounding residue survives — against scores of order
0.1 they are zero. Every meaningful digit of every other field matches, and the
other 21 outputs, including `.matched`, `.match_report` on both keys, the merged
matrix, `.pgs_corr`, `.best_pgs` and every MHC output, are byte-identical.

It is a property of the toolchain, not of any code change. The figures move a
little between releases because they depend on which commands the sweep runs,
not on the size of the error: an earlier, narrower measurement over the
partitioned `.sscore` alone reported 2 of 900 cells. Quote the sweep the numbers
came from, or they read as a regression when the sweep merely grew.

## Runtime requirements

Neither binary needs anything installed.

`adjutant-linux-x86_64-static` is verified by running it in `busybox`, which
carries no glibc whatsoever. That costs size — 28 MB against about 1.6 MB for a
dynamic build — and the bulk is OpenBLAS's per-CPU kernels, not debug info:
`strip --strip-unneeded` recovers only about 1 MB, so it is left unstripped.

`adjutant-macos-arm64` links only `Accelerate.framework`, `libSystem.B.dylib`,
`libz.1.dylib` and `libc++.1.dylib`, all shipped with the OS. Check with
`otool -L`; anything outside `/System` or `/usr/lib` means the build picked up a
Homebrew library and will not run on a machine without it.

The Linux binary is built on Debian bullseye (glibc 2.31, 2020) so it runs on
older long-lived cluster images.

## Notable changes to check against

**`--give-best-pgs cols=` pins the columns and their order.** `cols=PGS_ID,
PHENO,R2_INCR,P` emits exactly those, in that order, so a positional parser
stays correct across versions: adding a column changes the default layout,
never a named one. Prefer it over indexing the default layout by position.

**`.best_pgs` is now 20 columns, not 21.** `--folds` and its `FOLDS` column
have been removed: cross-validation corrected an in-sample optimism of about
1/(N-p) -- one parameter's worth, since `R2_INCR` is a difference of two R2 that
differ by one predictor -- and overshot it, landing further from the truth than
doing nothing at N=150. It did not address the winner's curse, which is the
bias that actually matters when ranking, and did not improve which model the
ranking picks. Anything parsing `.best_pgs` by column number needs updating;
by name, nothing changes except that `FOLDS` is gone.

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
- **`--sscore` no longer caps input width at 1024 columns.** The three row
  parsers tokenised into a fixed 1024-entry array, so a `.sscore` past 1024
  columns (about 1021 models, well within reach of a merged
  `--pgs-models-list` matrix) was rejected with "too many columns", and a
  `--pheno` past 1024 phenotypes or a `--covar` past 1023 covariates failed
  *silently* — exit 6 with nothing on stderr. The buffer is now sized from the
  header, once per file rather than once per sample.
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

**Purge `*.o` from the staged copy before building in the container.** A native
`make release` leaves 68 Mach-O objects in the tree. Copy those into the
container and `make` treats them as current, so the Linux link is handed macOS
objects. Excluding `bin/` is not enough — the objects sit beside the sources.

**A bind mount from an unshared path mounts empty and says nothing.** Docker
Desktop only shares configured prefixes; a scratchpad under `/private/tmp` is
not one. `-v /private/tmp/...:/src` produced an empty `/src`, the build did no
work, and the run exited 0 with no output at all. Stage under the home
directory, and check `ls /src | wc -l` inside the container before trusting a
build. Capture docker's own exit status into a variable — piping its output
through `tail` reports the exit code of `tail`.

**Cross-compiling to x86-64 on Apple Silicon is unreliable.** Under qemu,
`cc1`/`g++` take a SIGSEGV at random — a different translation unit each time,
on files that compile fine natively. `make` resumes from the object files
already produced, so retrying converges; drop to `-j1` after the first failure.
Verify the result from the ELF header (`7f454c46` magic, `e_machine` `0x3e`)
rather than with `file`, which is **not installed** in `debian:bullseye` — a
check that silently rejected 19 good builds. And exclude the host's `bin/`
from whatever you copy into the container: a stray macOS binary landing there
turns an existence check into a false pass.
