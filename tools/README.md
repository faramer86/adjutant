# tools/

Scaffolding used to reduce PLINK 2 down to `--score`. Kept because it is
reusable if the reduction is ever taken further; not needed to build or run
adjutant.

- `tools_strip.py` — brace-matched removal of if/else-if branches and whole
  function definitions from C++ sources. Refuses ambiguous matches and aborts
  if a deletion would unbalance braces.
- `reduce_plink2.py` — compiler-guided pass that strips `Plink2Cmdline` fields
  and the flag-parse branches belonging to removed commands.
- `autofix.py` — iteratively deletes code referencing removed commands, guided
  by compiler errors. Has a `PROTECTED` list of load-path symbols it must never
  touch; without it, it will happily delete `LoadPvar()`.

All three operate on the source tree in place. The golden harness
(`../score_check.sh`) is what proves they did no damage.

## Dead-code elimination

Five detectors, each for a category the others cannot see. Run them in a loop
until every one comes back empty — each pass strands work for the next.

- `deadcode.py` — unreachable top-level functions in a `.cc`. Blanks comments
  and literals first, then matches bare identifiers, so it sees function
  pointers without being fooled by prose. Header *inline bodies* count as
  roots; header *prototypes* do not.
- `closure.py` — reachability across a `.cc`/`.h` **pair**, from the symbols
  the rest of the tree actually names. Catches what per-file analysis cannot: a
  ring of dead functions that reference each other. `adjutant_matrix` was 33
  definitions of which 3 were reachable, and `deadcode.py` called it "0 dead".
- `deadheader.py` — unreferenced inline definitions in a `.h`.
- `deadtype.py` / `prune_types.sh` — unreferenced typedef'd structs, enums and
  flagsets. Removing a command strands its thread-context struct and its
  modifier flagset; those in turn *reference* other dead symbols, which is what
  keeps them looking alive to the function-level pass.
- `deaddecl.py` — prototypes whose definitions are gone. Deleting a function is
  two edits and the second gets missed; a declaration for a function that does
  not exist is legal right up until someone calls it, so these accumulate
  silently. `adjutant_matrix.h` had twelve.
- `linkmap.py` — asks the linker instead of guessing. Links twice, `-O2` and
  `-O0 -fno-inline`, with `-Wl,-dead_strip -Wl,-map`, and diffs the two
  dead-symbol reports. Dead in **both** means genuinely uncalled; dead at `-O2`
  **only** means either inlined (keep) or reachable solely through a branch the
  optimiser folded away (delete). The second bucket is a question, not an
  answer: `Count12Vec6` and `RandomizeBigstack` land in it together and only
  one is dead.
- `extroots.py` — roots for a module: symbols other objects actually reference,
  read from `nm`.
- `prune_named.py` — applies an explicit kill list, for when the verdict comes
  from `linkmap.py` or a human rather than from a reachability pass.
- `prune_module.sh` / `prune_include.sh` / `prune_header.sh` — strip, then
  **clean-build** + `score_check.sh` + `input_matrix.sh`, commit on green,
  revert on red.
- `verify_and_commit.sh` — the same gate for hand-edited changes. Refuses to
  run with uncommitted `tools/` changes: on red it does `git checkout -- .`,
  which will otherwise eat the fix for whatever bug the run just exposed.
- `test_cli_guards.py` — asserts that a mistyped mode never writes to a source
  file, and that a bare invocation prints usage instead of a traceback. Run it
  after touching any argument handling: `python3 tools/test_cli_guards.py`.

The clean build is not optional. An incremental build leaves stale `.o` files
for TUs that include a pruned header, so the link succeeds on symbols that no
longer exist and the whole check passes against a binary that cannot be rebuilt.

`CORE.md` records the protected call chain these tools must never break.

### Ways these tools have been wrong

Every one of these shipped a broken build before the harness caught it. They
are the reason the clean-build gate exists; do not skip it.

- **A dead call site still counts as a reference.** `RandomizeBigstack` is
  called under `if (randmem)`, where `randmem` is a literal `0` that is never
  assigned. Greps cannot see this; the linker can.
- **Macro bodies are call sites.** `BigstackAllocX` is reached only through
  `#define BIGSTACK_ALLOC_X(...)`. And the directive may be written `#  define`
  — match the directive, not the literal prefix.
- **Index every `#ifdef` arm.** `plink2_bits.cc` defines
  `PopcountWordsIntersect` twice, once per arm. Keeping only the first meant
  walking the AVX2 body and never the one that actually compiles, so its callee
  looked dead. This inflated `plink2_simd` to "68 dead" when the true figure
  was 2.
- **Anonymous typedefs have no tag.** `deadtype.py` assumed a name appears
  twice in its own definition; `typedef union { ... } CmpExprArgs;` has it once,
  so the threshold silently covered one real use. Excise the span and search
  what remains instead of counting.
- **A type can be dead while its members are not.** `typedef enum { kChrOffsetX,
  ... } ChrOffset;` — the type name appears nowhere else, the constants are
  everywhere. Same for `FLAGSET_DEF`, whose members are named `kf<Name>*`.
- **A mistyped mode used to mean "strip".** Every mode-taking script put the
  destructive path in the *fallback* — `if mode == 'report': ... return` with
  the strip code below it, an unguarded `else:`, or `cut_branch if mode ==
  'flag' else cut_function`. So any string that was not the read-only mode
  wrote to the file. `deadcode.py list adjutant_pgs_eval.cc ""` printed
  "removed 5/5 dead functions" and took a 3054-line file down to 1280, with no
  backup and no confirmation; recovery was by hand. Modes are now validated
  against the documented set up front, as `linkmap.py` always did. Keep the
  destructive path behind an explicit mode, never behind a fallback.
- **A missing argument used to mean "go".** The two scripts with no mode to
  mistype were left out of that fix and had the same hole one level down:
  `autofix.py` defaulted its target (`sys.argv[1] if len(sys.argv) > 1 else
  'adjutant.cc'`) and `reduce_plink2.py` read no arguments at all and hardcoded
  the name. So a bare `python3 tools/autofix.py` — a stray shell line, a
  half-typed command, a tab-completion that fired early — put `adjutant.cc`
  through 30 rounds of compiler-guided deletion and wrote the result back.
  Both now want the invocation spelled out: a path for `autofix.py`, `--apply`
  for `reduce_plink2.py`, usage and exit 1 otherwise. Destructive-by-default is
  the same bug as destructive-by-fallback; the argument that selects the write
  has to be one the user actually typed.
- **Always keep `deadcode.py`'s safety net.** `index_functions` does not
  recognise every signature style, so a caller can go unindexed and its body
  never scanned. Before deleting anything, confirm the name really does not
  occur outside its own definition. Adding this net to `closure.py` took
  `plink2_bits` from "12 dead" to 0 — all twelve were false positives.

## Distribution builds

Not part of the four suites. These produce the binaries under `dist/`, which is
gitignored; two targets are supported, macOS arm64 and Linux x86-64.

- `build_linux.sh <docker-platform> <output-name>` — dynamic build inside a
  Debian bullseye container. glibc 2.31 (2020) is the baseline on purpose, so
  the result runs on older long-lived cluster images. `libstdc++`/`libgcc` are
  linked statically; OpenBLAS stays dynamic, so the target needs `libopenblas0`
  (which pulls in `libgfortran5`).
- `build_linux_static.sh <docker-platform> <output-name>` — the same, fully
  static: no glibc, no OpenBLAS, nothing. Costs ~27 MB against ~1.5 MB, which is
  OpenBLAS's per-CPU kernels rather than debug info.
- `sscore_cmp.sh <golden-dir> <out-dir>` — compares `.sscore` files and, where
  they differ, reports **how big** the differing numbers are. "Not byte-identical"
  and "numerically wrong" are different claims: a value of 1e-18 where the truth
  is 0 is cancellation residue, and only the magnitude separates the two.
  Both build scripts call it to report their own verification.

macOS is just `make release`; there is no script for it.

### Four things these encode that are easy to get wrong

- **`rm -f` before `cp`.** `cp` overwrites in place, and overwriting the text
  pages of a running executable wedges that process in uninterruptible state on
  macOS — after which every later exec of the same inode wedges too, at zero CPU,
  looking exactly like an application hang. The wedged processes cannot be killed.
- **Static link order.** With `-static` the linker resolves left to right:
  `-llapacke -lopenblas -lgfortran [-lquadmath]`. `libquadmath` exists on x86-64
  but not aarch64.
- **colima shares only paths under `$HOME`.** A `/tmp` bind mount arrives
  silently empty, so `sscore_cmp.sh` is passed through the environment as base64
  rather than mounted, and output goes to `<repo>/dist`.
- **`score_check.sh` cannot pass unmodified on Linux.** `score_regression.sh`
  normalises only the volatile fields macOS emits; Linux adds `, ~NNNN available`
  to the RAM line and a `Using up to N compute threads.` line. That is log
  formatting, not scoring, so the `.sscore` comparison carries the verdict.

## Differential testing against stock plink2

- `score_diff_vs_plink2.sh` — runs ~50 `--score` invocations through both
  adjutant and a stock `plink2` on data generated by `plink2 --dummy` with fixed
  seeds, and compares the `.sscore` bytes. Needs `plink2` on PATH; writes its
  filesets to `/tmp/adjutant_score_diff` (override with `WORK_DIR`), never into
  the repo.

Unlike `score_check.sh`, which compares against committed goldens, this checks
against plink2 *itself* on data neither has seen. That is what found the open
`se` discrepancy: under the `se` modifier, variants at extreme MAF — the ones
that engage the sparse/difflist genotype path — disagree between adjutant and
plink2 a.6.9. The golden dataset has no variants below MAF 0.01, so
`score_check.sh` cannot see it. Attribution needs a pristine a.7.1.a binary; see
the memory note on external oracles.
