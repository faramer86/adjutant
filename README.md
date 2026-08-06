# Adjutant 🫡

![status: beta](https://img.shields.io/badge/status-beta-orange)
![license: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue)
![platform: macOS | Linux](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-lightgrey)
![derived from PLINK 2](https://img.shields.io/badge/derived%20from-PLINK%202.0-informational)

> [!WARNING]
> **Beta.** Adjutant is under active development and issues may arise. The
> `--score` engine is unmodified upstream code, checked against stock plink2 on
> every build, but the flags built on top of it are newer and less exercised in
> the wild. Check your results, and please
> [open an issue](https://github.com/faramer86/adjutant/issues) if something
> looks wrong — the command line and the log are enough to act on.

Adjutant is **PLINK 2** reduced to the functions for polygenic-score
experiments, and those functions extended. It aligns scores to the provided
reference, it computes them, it evaluates them (e.g. Incremental R2) — you only
need to give an order.

One binary, no runtime dependencies, quick, easy, with `--score` output
byte-identical to stock plink2.

> **Most of this code is reduced [PLINK 2.0](https://github.com/chrchang/plink-ng)'s**,
> and the credit belongs there. `adjutant` was made by deleting commands and reusing functions to create new functions, not by writing a scorer. If you find it useful, hit a star for this repository and cite PLINK 2 — see [NOTICE](NOTICE).

---

## Three modes

| | Input | Output | Flags |
|---|---|---|---|
| **Score** | genotypes + weights | `.sscore` | `--score` |
| **Evaluate** | scores + a phenotype | `.best_pgs`, `.pgs_corr` | `--give-best-pgs`, `--give-corr-matrix` |
| **Harmonise** | weights + a `.bim` | `.matched` | `--pgs-model`, `--pgs-models-list` |

Score and evaluate can run in one command; evaluate also takes a score matrix
you already have (`--sscore`). Harmonise is standalone — two text files in, two
out, no genotypes.

---

## Build

```sh
make release          # -> bin/adjutant
```

Needs a C++17 compiler, zlib, and a BLAS/LAPACK (Accelerate on macOS, OpenBLAS
elsewhere). Everything else is vendored.

Or take a pre-built binary. Neither needs anything installed:

| platform | download |
|---|---|
| Linux x86-64 | [`adjutant-linux-x86_64-static`](https://github.com/faramer86/adjutant/releases/latest/download/adjutant-linux-x86_64-static) |
| macOS Apple Silicon | [`adjutant-macos-arm64`](https://github.com/faramer86/adjutant/releases/latest/download/adjutant-macos-arm64) |

```sh
curl -LO https://github.com/faramer86/adjutant/releases/latest/download/adjutant-linux-x86_64-static
chmod +x adjutant-linux-x86_64-static
./adjutant-linux-x86_64-static --version
```

The Linux build is fully static — it runs on anything x86-64, including images
with no glibc. The macOS build links only frameworks that ship with the OS.
Checksums are in [`SHA256SUMS`](https://github.com/faramer86/adjutant/releases/latest/download/SHA256SUMS)
on the same release.

---

## Quick start

**Every example below runs against [`example/`](example/)** — a complete,
committed set of working input files. Clone the repo and they run as written:

```sh
cd example
```

300 samples and 40 variants, small enough to read, and not random noise:
`HEIGHT` is generated from `PGS_B`, so the tool has a real answer to find. The
genotypes ship as both a `.bed/.bim/.fam` trio and a `.pgen/.pvar/.psam` one;
the examples use `--bfile`, and `--pfile demo` works identically.

### The input files, and what each is for

| File | Needed by | Format |
|---|---|---|
| `demo.bed/.bim/.fam` | `--bfile` | genotypes, PLINK 1 binary |
| `demo.pgen/.pvar/.psam` | `--pfile` | genotypes, PLINK 2 |
| `weights.txt` | `--score` | variant ID, effect allele, then one weight column per model |
| `pheno.txt` | `--pheno` | `#IID` then one column per trait |
| `covar.txt` | `--covar` | `#IID` then one column per covariate |
| `partitions.txt` | `--score-partitions` | variant ID, then one 0/1 column per partition |
| `pairs.txt` | `--pgs-pheno-pairs` | two columns: phenotype, PGS |
| `target.bim` | `--bim` | the target's allele coding |
| `published_pgs.txt` | `--pgs-model` | chr, pos, effect allele, other allele, beta, freq |
| `demo.afreq` | `--read-freq` | allele frequencies |

Copy any as a template for your own data — see
[example/README.md](example/README.md) for how they were made.

---

## I. Scoring

### One model

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --out result
```

The numbers are **column positions in the weight file**: variant ID, effect
allele, weight. Each defaults to the one after the last, so `1 2` means `1 2 3`
— the weight column is chosen for you.

```
--score <file> <id> <effect allele> <weight>
```

Any of them can sit anywhere — name all three when your file is laid out
differently. `weights.txt` holds `PGS_C` in column 5, so this scores that model
instead:

```sh
adjutant --bfile demo --score weights.txt 1 2 5 header-read --out result
```

`header-read` takes score names from the first line.

### Several models in one pass

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read \
         --score-col-nums 3-5 --out result
```

`result.sscore`:

| FID | IID | ALLELE_CT | NAMED_ALLELE_DOSAGE_SUM | PGS_A_AVG | PGS_B_AVG | PGS_C_AVG |
|---|---|---|---|---|---|---|
| 0 | per0 | 80 | 34 | -0.0175 | 0.01375 | 0.00125 |
| 0 | per1 | 80 | 34 | -0.03 | 0.0375 | -0.01 |
| 0 | per2 | 80 | 48 | -0.00375 | 0.02125 | -0.01375 |

The `FID` column comes from the `.fam`, which always carries one; `--pfile`
output on a `.psam` without FIDs has no such column. Everything else is
identical.

### Split a model into partitions

`--score-partitions` takes a 0/1 mask per variant and scores each subset
separately in the same pass. A partition is any set of variants you want scored
on its own; the file defines them, the header names them:

```
SNP    PART_A  PART_B
snp0   1       1
snp1   0       1
```

A variant may belong to several partitions, or to none. Values must be exactly
`0` or `1`.

```sh
adjutant --bfile demo --score weights.txt 1 2 4 header-read \
         --score-partitions partitions.txt --out result
```

| FID | IID | ALLELE_CT | NAMED_ALLELE_DOSAGE_SUM | PGS_B_AVG | PGS_B_PART_A_AVG | PGS_B_PART_B_AVG |
|---|---|---|---|---|---|---|
| 0 | per0 | 80 | 34 | 0.01375 | 0.0075 | 0.005 |
| 0 | per1 | 80 | 34 | 0.0375 | 0.04 | 0.0525 |
| 0 | per2 | 80 | 48 | 0.02125 | -0.0025 | 0.045 |

`PGS_B_AVG` is the unpartitioned total. Each partition divides by **its own**
variants' allele count, so partitions of different sizes stay comparable — and
therefore do **not** sum to the total. `cols=+scoresums` gives sums, which do.

### Drop the MHC

Long-range LD and extreme polymorphism make the MHC unlike the rest of the
genome, so PGS work often excludes it. Pick the flag for **your data's build** —
nothing auto-detects, and the wrong flag silently drops the wrong window:

| flag | region | chr6 interval |
|---|---|---|
| `--exclude-mhc-hg19` | MHC | 28,477,797 – 33,448,354 |
| `--exclude-mhc-hg38` | MHC | 28,510,120 – 33,480,577 |
| `--exclude-xmhc-hg19` | extended MHC | 25,726,291 – 33,368,421 |
| `--exclude-xmhc-hg38` | extended MHC | 25,726,063 – 33,400,644 |

Bounds are inclusive. At most one flag per run.

```sh
adjutant --bfile mhc --exclude-mhc-hg19 --score mhc_weights.txt 1 2 3 --out result
```

```
--exclude-mhc-hg19: chr6:28477797-33448354, 10 variants excluded.
--score: 14 variants processed.
```

Every run logs the interval it applied. In harmonisation the flags apply to
**both** the `.bim` and the weight file, so excluded variants are absent from
the `.match_report` rather than reported as unmatched:

```sh
adjutant --bim mhc.bim --pgs-model mhc_pgs.txt 1 2 3 4 5 --exclude-xmhc-hg38 --out h
```

> [!IMPORTANT]
> **The extended region does not contain the plain one.** They come from
> different published definitions: the GRC MHC region ends **79,933 bp past** the
> xMHC end, in both builds. So `--exclude-xmhc-hg19` reaches 2.75 Mb further 5′
> but stops short 3′, and can exclude **fewer** variants than
> `--exclude-mhc-hg19` — on the `example/mhc` fileset, 9 against 10.
>
> No published MHC/xMHC pairing nests, so this is intrinsic to the definitions,
> not a quirk of this tool: adjutant reports each region exactly as published
> rather than silently widening one to fit the other. To exclude everything, run
> the plain flag in a separate pass and intersect.

---

## II. Evaluation

### Rank models against a trait

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
         --pheno pheno.txt --pheno-name HEIGHT --covar covar.txt \
         --give-best-pgs --out result
```

`result.best_pgs`, sorted by `R2_INCR`:

| PGS_ID | PHENO | TYPE | N | BETA | SE | P | R2_INCR | R2_INCR_L95 | R2_INCR_U95 |
|---|---|---|---|---|---|---|---|---|---|
| PGS_B | HEIGHT | QT | 300 | 1.04992 | 0.133779 | 7.77e-14 | 0.171924 | 0.099872 | 0.254388 |
| PGS_C | HEIGHT | QT | 300 | 0.0495218 | 0.147079 | 0.736581 | 0.000382137 | 0 | 0.0176273 |
| PGS_A | HEIGHT | QT | 300 | 0.00247 | 0.147029 | 0.986608 | 9.51667e-07 | 0 | 0.0131051 |

PGS_B is the model `HEIGHT` was built from, and it wins. The two noise models
get an `R2_INCR_L95` of exactly **0** — their interval contains zero. Read the
interval before the ranking.

`BETA` is per standard deviation of the score. Two models are fitted per row:
phenotype on covariates alone, then with the score added.

### Case/control traits

Detected automatically. Add `--prevalence` for liability-scale R²:

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
         --pheno pheno.txt --pheno-name DISEASE --covar covar.txt \
         --give-best-pgs --prevalence 0.05 --out result
```

| PGS_ID | PHENO | TYPE | N | N_CASE | OR | P | R2_INCR | R2L_INCR | AUC_NULL | AUC_FULL |
|---|---|---|---|---|---|---|---|---|---|---|
| PGS_B | DISEASE | CC | 300 | 148 | 2.38262 | 7.76e-09 | 0.128126 | 0.108728 | 0.541474 | 0.701814 |
| PGS_A | DISEASE | CC | 300 | 148 | 0.920839 | 0.478059 | 0.00167475 | 0.0014212 | 0.541474 | 0.540985 |
| PGS_C | DISEASE | CC | 300 | 148 | 1.04041 | 0.732916 | 0.000387154 | 0.000328542 | 0.541474 | 0.54183 |

### Cross-validate the ranking

Picking the best of N models and reporting its in-sample R² overstates it.
`--folds k` fits on k−1 folds and scores on the held-out one:

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
         --pheno pheno.txt --pheno-name HEIGHT --covar covar.txt \
         --give-best-pgs --folds 5 --out result
```

| PGS_ID | FOLDS | R2_NULL | R2_FULL | R2_INCR | R2_INCR_L95 |
|---|---|---|---|---|---|
| PGS_B | 5 | -0.0228002 | 0.14881 | 0.171611 | 0.0985724 |
| PGS_C | 5 | -0.0228002 | -0.0268428 | -0.00404261 | NA |
| PGS_A | 5 | -0.0228002 | -0.0298362 | -0.00703593 | NA |

The real signal survives at 0.172. The two noise models go **negative** — they
predict worse than the covariates alone — and a negative R² gets `NA` bounds.
That is the flag earning its keep.

`--folds 1` is the default, bit-identical to omitting the flag.

### Correlate the models

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
         --give-corr-matrix --out result
```

| ID | PGS_A | PGS_B | PGS_C |
|---|---|---|---|
| PGS_A | 1 | -0.0248890537084 | -0.0589350562717 |
| PGS_B | -0.0248890537084 | 1 | -0.0752174055941 |
| PGS_C | -0.0589350562717 | -0.0752174055941 | 1 |

Needs no phenotype. Read it with
`read.table(f, header=TRUE, comment.char="", row.names=1)` in R, or
`read_csv(f, sep='\t', index_col=0)` in pandas.

### Evaluate only certain pairs

When each score was built for one specific trait, the full cross-product is
mostly noise. `pairs.txt`:

```
# phenotype  PGS
HEIGHT       PGS_B
DISEASE      PGS_B
```

```sh
adjutant --sscore result.sscore --pheno pheno.txt --covar covar.txt \
         --give-best-pgs --pgs-pheno-pairs pairs.txt --out result
```

| PGS_ID | PHENO | N | R2_INCR | P |
|---|---|---|---|---|
| PGS_B | HEIGHT | 300 | 0.171924 | 7.77e-14 |
| PGS_B | DISEASE | 300 | 0.128126 | 7.76e-09 |

Two rows instead of six. An unknown name is an error, not a silently dropped
row.

### A score matrix you already have

`--sscore` takes any `.sscore` from plink2 or adjutant, and needs no genotypes:

```sh
adjutant --sscore result.sscore --pheno pheno.txt --covar covar.txt \
         --give-best-pgs --out ranked
```

---

## III. Harmonisation

A published PGS rarely uses your dataset's allele coding. `--pgs-model` aligns
it, reimplementing `bigsnpr::snp_match()`.

```sh
adjutant --bim target.bim --pgs-model published_pgs.txt 1 2 3 4 5 --out harmonised
```

The five numbers are the columns holding chromosome, position, effect allele,
other allele and weight.

```
--pgs-model: 40 source variants, 40 found in the .bim by position.
--pgs-model: 14 ambiguous (A/T, C/G) dropped; use --allow-ambiguous to keep them.
--pgs-model: 26 matched (8 strand-flipped, 10 reversed), 0 unmatched on alleles.
```

`harmonised.matched` carries the **target's** identity — chromosome, position,
both alleles — taken from the matched `.bim` row, so the file stands on its own:

| CHR | POS | ID | A1 | A2 | BETA |
|---|---|---|---|---|---|
| 1 | 2 | snp2 | A | G | 0.4 |
| 1 | 3 | snp3 | T | C | 0.4 |
| 1 | 4 | snp4 | C | T | -0.1 |

It feeds straight back into `--score`, with `ID` in column 3 and `A1` in 4:

```sh
adjutant --bfile demo --score harmonised.matched 3 4 header-read --score-col-nums 6 --out result
```

`harmonised.match_report` keeps every variant, including the dropped ones:

| STATUS | count | meaning |
|---|---|---|
| `EXACT` | 8 | alleles already match |
| `REVERSED` | 10 | effect/other swapped — **beta negated** |
| `FLIPPED` | 8 | opposite strand |
| `AMBIGUOUS` | 14 | A/T or C/G — strand undecidable, dropped |

**Weights must be betas.** Negating an odds ratio is meaningless; adjutant warns
if every weight is positive with mean in [0.5, 2.0].

### Matching on rsID instead of position

If the weight file and your target come from **different genome builds**,
`chr:pos` matching finds nothing. A leading `id` joins on the variant ID
instead:

```sh
adjutant --bim target.bim --pgs-model pgs_rsid.txt id 1 4 5 6 --out harmonised
```

Four numbers now, not five — ID, effect allele, other allele, weight — since
chromosome and position go unused. `example/pgs_rsid.txt` is exactly this case:
an rsID-keyed PGS whose positions belong to another build.

```
--pgs-model: 40 source variants, 40 found in the .bim by ID.
--pgs-model: 100.0% of source variants found in the .bim by ID.
```

Matched on position, the same file finds **nothing at all**:

```
--pgs-model: 40 source variants, 0 found in the .bim by position.
```

Everything else is unchanged: allele matching, beta negation, the report and
`--exclude-mhc-*` behave identically, and `.matched` still carries the target's
coordinates.

> [!IMPORTANT]
> **This trades a build problem for an rsID-versioning one.** Merged and retired
> rsIDs mean a source and target built against different dbSNP releases can
> silently fail to match, or match the *wrong* variant — and unlike a build
> mismatch, invisible in the coordinates, which never enter the comparison.
>
> So the run prints the match rate: a join that found 12% means the two files
> disagree about what an rsID means. Read it before you trust the output.
>
> Rows with a `.` ID are skipped on both sides. Where the `.bim` repeats an ID
> the first row wins; unreachable rows are listed in `<prefix>.duplicate_ids`.
> Consider re-IDing both files as `chr_pos_ref_alt`, unique by construction.

### Rescuing strand-ambiguous variants

A/T and C/G variants look the same on either strand. Given frequencies from a
**comparable ancestry** on both sides, `--resolve-ambiguous` decides them:

```sh
adjutant --bim target.bim --pgs-model published_pgs.txt 1 2 3 4 5 6 \
         --read-freq demo.afreq --resolve-ambiguous --out resolved
```

```
--pgs-model: 14 ambiguous; 11 resolved by frequency (--resolve-ambiguous, MAF < 0.4 on both sides),
             3 unresolved and dropped.
--pgs-model: 37 matched (8 strand-flipped, 13 reversed), 0 unmatched on alleles.
```

26 matched becomes 37. The 3 left over have a MAF too near 0.5 to call and are
refused rather than guessed. Column 6 is the source effect-allele frequency.

### Merging several models

`--pgs-models-list` harmonises several published models against the **same**
target and merges them into one matrix — one column per model, one row per
variant matched by at least one of them.

Give it a text file listing the weight files, one per line:

```
pgs_alpha.txt
pgs_beta.txt
```

Blank lines and `#` comments are skipped; paths are relative to where you run
the command. Each column takes its name from the file's basename up to the first
dot, so `PGS000123.txt.gz` becomes `PGS000123`; two paths yielding the same name
is an error.

Every model must use the **same column layout**; the five indices apply to all
of them:

```
chr  pos  effect_allele  other_allele  beta
1    14   T              C             -0.50
1    15   G              T             0.20
```

```sh
adjutant --bim target.bim --pgs-models-list models.txt 1 2 3 4 5 --out merged
```

```
--pgs-models-list: pgs_alpha: 25 source variants, 16 matched (0 strand-flipped, 0 reversed).
--pgs-models-list: pgs_beta:  26 source variants, 15 matched (0 strand-flipped, 5 reversed).
--pgs-models-list: 2 models, 26 variants in the union.
```

`merged.matched` is directly usable as a `--score` input:

| ID | A1 | pgs_alpha | pgs_beta |
|---|---|---|---|
| snp2 | A | 0.4 | **0** |
| snp14 | C | -0.4 | 0.5 |
| snp19 | T | **0** | 0.4 |

A model that does not carry a variant gets `0`, not a missing value. Of the 26
variants in the union only 3 carry a nonzero weight from *both* models: the
files overlap on 11 source rows, but most are strand-ambiguous (A/T, C/G) and
dropped unless you pass `--allow-ambiguous` or `--resolve-ambiguous`. So a zero
means the model lacked the variant *or* its match was rejected —
`.match_report` tells them apart.

**Harmonising before merging is the whole point.** Two models can name opposite
effect alleles for the same variant; merging raw weight files on variant ID
would put one column's beta on the wrong allele, and nothing downstream would
notice. `merged.match_report` shows the reorientation:

| #MODEL | CHR | POS | SRC_EFFECT | SRC_OTHER | SRC_BETA | STATUS | BIM_ID | OUT_A1 | OUT_BETA |
|---|---|---|---|---|---|---|---|---|---|
| pgs_beta | 1 | 14 | T | C | -0.5 | **REVERSED** | snp14 | C | **0.5** |

`pgs_beta` measured its effect on `T`, the target's A1 is `C`, so the weight is
re-expressed as `+0.5` on `C`. Every column in `.matched` is against the same
A1, which is what makes them comparable.

> [!NOTE]
> **Reorienting a weight shifts absolute scores by a constant.** For a biallelic
> variant `dosage(other) = 2 − dosage(effect)`, so
>
> `β · dosage(effect) = 2β − β · dosage(other)`
>
> Writing `−β` against the other allele reproduces the original contribution
> **minus `2β`** — over the reversed variants, one constant identical for every
> sample. Anything fitted with an intercept absorbs it, so **rankings,
> correlations, R² and regression coefficients are unchanged**; only raw
> `SCORE_SUM` moves. On 5 variants at |β| = 0.4, every sample's sum differed by
> exactly `4.0` (= 2 × 0.4 × 5), `r = 1.000000000000` across 200 samples.
>
> Exact for **additive** models only. Dominant and recessive effect types break
> the identity, so tools supporting them must split a model set across files
> instead of reorienting it; adjutant is additive-only, so one merged file
> always suffices.

---

## Output reference

### `.sscore`

| Column | Meaning |
|---|---|
| `IID` (and `FID`) | sample identifiers |
| `ALLELE_CT` | non-missing allele count behind the score |
| `NAMED_ALLELE_DOSAGE_SUM` | summed dosage of the named alleles |
| `<model>_AVG` | score, averaged over alleles |
| `<model>_SUM` | score, unaveraged (`cols=+scoresums`) |

### `.best_pgs`

21 columns, the same set for both trait types.

| Column | Meaning |
|---|---|
| `PGS_ID`, `PHENO` | which model, which phenotype |
| `TYPE` | `QT` or `CC` |
| `N`, `N_CASE` | complete cases; cases among them (`CC` only) |
| `FOLDS` | cross-validation folds behind `R2_INCR` |
| `BETA`, `SE`, `P` | per SD of the score. `CC`: log-odds from a logistic fit |
| `OR` | `exp(BETA)` (`CC` only) |
| `R2_NULL`, `R2_FULL` | variance explained without and with the score |
| `R2_INCR` | their difference — **the sort key** |
| `R2_INCR_L95` / `_U95` | 95% CI, Fisher z on the partial correlation |
| `R2L_INCR` + CI | liability scale (`CC`, needs `--prevalence`) |
| `AUC_NULL`, `AUC_FULL`, `AUC_INCR` | `CC` only; always in-sample |

`R2_INCR` is **observed-scale for both trait types** — for case/control it comes
from an OLS fit on the 0/1 outcome, not a pseudo-R² of the logistic fit, because
that is what the liability transform assumes.

Non-converged rows sort last, `NA` in every numeric field.

---

## Checking a build

`--score` output is **byte-identical to stock plink2**. That is the design
contract, and why the build uses `-ffp-contract=off`: the flag stops the compiler
fusing multiply-add in the scoring dot product, which keeps `-O3` identical to
`-O2` and to upstream.

You need not take that on trust. If you have plink2 installed, check against
your own data:

```sh
plink2   --bfile mydata --score w.txt 1 2 header-read --out ref
adjutant --bfile mydata --score w.txt 1 2 header-read --out mine
diff ref.sscore mine.sscore        # expected: no output
```

---

## Provenance and licence

Derived from `PLINK v2.0.0-a.7.1.a`. The scoring engine, the `.pgen`/`.pvar`/
`.psam` readers, the filters, and the text and compression layers are upstream
code, modified little or not at all.

- **[LICENSE](LICENSE)** — GPL v3, the application.
- **[LICENSE.LESSER](LICENSE.LESSER)** — LGPL v3, the libraries under `include/`
  that upstream ships separately.
- **[NOTICE](NOTICE)** — full attribution and the vendored third-party licences.

Cite PLINK 2:

> Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015)
> Second-generation PLINK: rising to the challenge of larger and richer
> datasets. *GigaScience* 4.
