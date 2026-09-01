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
> looks wrong - the command line and the log are enough to act on.

**Adjutant** began with a simple idea: reduce the **PLINK 2** codebase to the core functionality required for polygenic score (PGS) analysis, 
then extend it with additional functions commonly needed in PGS experiments. 
It can calculate, evaluate, harmonise, merge polygenic scores through a familiar plink–style command-line interface.
It is distributed as a single binary with no runtime dependencies.
It is quick, easy, with `--score` output byte-identical to stock PLINK 2.

> **Most of this code is reduced [PLINK 2.0](https://github.com/chrchang/plink-ng)'s**,
> and the credit belongs there. `adjutant` was made by deleting commands and reusing functions to create new functions, not by writing a scorer. 
> If you find it useful, hit a star for this repository and cite PLINK 2 — see [NOTICE](NOTICE).

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
---

## Quick start

**Every example below runs against [`example/`](example/)** — a complete,
committed set of working input files. Clone the repo and they run as written:

```sh
cd example
```

### The input files, and what each is for

| File | Needed by | Format |
|---|---|---|
| `demo.bed/.bim/.fam` | `--bfile` | genotypes, PLINK 1 input files |
| `demo.pgen/.pvar/.psam` | `--pfile` | genotypes, PLINK 2 input files |
| `weights.txt` | `--score` | variant ID, effect allele, then one weight column per model |
| `pheno.txt` | `--pheno` | `#IID` then one column per trait |
| `covar.txt` | `--covar` | `#IID` then one column per covariate |
| `partitions.txt` | `--score-partitions` | variant ID, then one 0/1 column per partition |
| `pairs.txt` | `--pgs-pheno-pairs` | two columns: phenotype, PGS |
| `target.bim` | `--bim` | the target's allele coding |
| `published_pgs.txt` | `--pgs-model` | chr, pos, effect allele, other allele, beta, freq |
| `demo.afreq` | `--read-freq` | allele frequencies |

See [example/README.md](example/README.md) for how they were made.

---

## What it can do

**[I. Scoring](#i-scoring)** — genotypes + weights → `.sscore`

- **[One model](#1-one-model)** — Apply a single weight file to a genotype fileset.
- **[Several at once](#2-several-models-in-one-pass)** — Score any number of weight columns in a single pass, one output column per model (`--score-col-nums`).
- **[Variant subsets](#3-split-a-model-into-partitions)** — Score chosen subsets of variants separately, alongside the unpartitioned total (`--score-partitions`).
- **[Drop the MHC](#4-drop-the-mhc)** — Exclude the chromosome 6 MHC region, using coordinates matched to your genome build (`--exclude-mhc-*`).

**[II. Evaluation](#ii-evaluation)** — scores + a phenotype → `.best_pgs`, `.pgs_corr`

- **[Rank models](#1-rank-models-against-a-trait)** — Rank candidate scores by the R² each adds to a covariates-only null model (`--give-best-pgs`).
- **[Case/control traits](#2-casecontrol-traits)** — Binary phenotypes are detected automatically and reported with odds ratios and AUC; `--prevalence` adds liability-scale R².
- **[Correlate models](#3-correlation-matrix)** — Compute a Pearson correlation matrix across the scored models (`--give-corr-matrix`).
- **[Chosen pairs only](#4-evaluate-only-certain-pairs)** — Restrict evaluation to named (phenotype, model) pairs instead of the full cross-product (`--pgs-pheno-pairs`).
- **[One PGS, every trait](#5-evaluate-single-pgs-against-all-of-the-phenotypes)** — Evaluate a score against every phenotype in the file by omitting `--pheno-name`.
- **[Reuse a score matrix](#6-use-a-score-matrix-you-already-have)** — Evaluate scores computed earlier, without reading genotypes again (`--sscore`).
- **[Pin the columns](#8-pin-the-columns-and-their-positions)** — Name the `.best_pgs` columns and their order, so a position-based parser stays correct across versions (`--give-best-pgs cols=`).
- **[Trim empty columns](#7-drop-the-columns-your-run-cannot-fill)** — Omit the `.best_pgs` columns a run cannot fill, which is 8 of 20 when every phenotype is quantitative (`--give-best-pgs no-na`).

**[III. Harmonisation](#iii-harmonisation)** — weights + a `.bim` → `.matched`

- **[Align a published PGS](#1-align-a-published-pgs-to-your-target-bim)** — Match a weight file to your target by position, flipping strand and negating effect sizes where required (`--pgs-model`).
- **[Match on rsID](#2-matching-on-rsid-instead-of-position)** — Join on variant ID rather than position, so a difference in genome build does not prevent matching (the `id` keyword).
- **[Rescue palindromes](#3-rescuing-strand-ambiguous-variants)** — Resolve strand-ambiguous A/T and C/G variants by comparing allele frequencies (`--resolve-ambiguous`).
- **[Merge many models](#4-merging-several-models)** — Harmonise several models against one target and merge them into a single matrix (`--pgs-models-list`).

---

## I. Scoring

### 1) One model

```sh
adjutant --bfile demo --score weights.txt 1 2 5 header-read --out result
```

The numbers are **column positions in the weight file**: variant ID, effect
allele, weight.

```
--score <file> <id> <effect allele> <weight>
```

`header-read` takes score names from the first line.

### 2) Several models in one pass

If you have more than one PGS to score use `--score-col-nums`.

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

### 3) Split a model into partitions

`--score-partitions` takes a 0/1 mask per variant and scores each subset
separately in the same pass. A partition is any set of variants you want scored
on its own; the file defines them, the header names them:

```
SNP    PART_A  PART_B
snp0   1       1
snp1   0       1
```

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
therefore do **not** sum to the total. `--score ... cols=+scoresums` gives sums, which do.

### 4) Drop the MHC

Long-range LD and extreme polymorphism make the major histocompatibility complex (MHC) locus unlike the rest of the
genome, so PGS work often excludes it. Pick the flag for **your data's build** -
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

It can be used for harmonisation step. 
In this case flags apply to **both** the `.bim` and the weight file, so excluded variants are absent from
the `.match_report` rather than reported as unmatched:

```sh
adjutant --bim mhc.bim --pgs-model mhc_pgs.txt 1 2 3 4 5 --exclude-xmhc-hg38 --out h
```

---

## II. Evaluation

### 1) Rank models against a trait

You can provide several scores, a phenotype and covariates, and `--give-best-pgs` ranks them by the variance explained by PGS.

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
         --pheno pheno.txt --pheno-name HEIGHT --covar covar.txt \
         --give-best-pgs --out result
```

`result.best_pgs`, sorted by `R2_INCR`:

| PGS_ID | PHENO | TYPE | N | BETA | SE | P | R2_NULL | R2_FULL | R2_INCR | R2_INCR_L95 | R2_INCR_U95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PGS_B | HEIGHT | QT | 300 | 1.04992 | 0.133779 | 7.77127e-14 | 0.00186797 | 0.173792 | 0.171924 | 0.099872 | 0.254388 |
| PGS_C | HEIGHT | QT | 300 | 0.0495218 | 0.147079 | 0.736581 | 0.00186797 | 0.00225011 | 0.000382137 | 0 | 0.0176273 |
| PGS_A | HEIGHT | QT | 300 | 0.00247 | 0.147029 | 0.986608 | 0.00186797 | 0.00186892 | 9.51667e-07 | 0 | 0.0131051 |

`BETA` is per standard deviation of the score. Two models are fitted per row:
phenotype on covariates alone (`R2_NULL`), then with the score added
(`R2_FULL`). `R2_INCR` is the difference, which is what the ranking sorts on —
`0.173792 − 0.00186797 = 0.171924`. `R2_NULL` is identical down the column
because the null model does not involve the score; it moves only when the
phenotype or the covariates change.

**The winner's `R2_INCR` is not an unbiased estimate of its performance.** The
ranking picks the largest of N noisy numbers, so the top row is upward-biased by
however much luck contributed — the smaller the gap between the leaders relative
to the standard error, the larger that bias. It is a property of taking a
maximum, not of the fit, so `R2_INCR_L95`/`_U95` do not correct for it either:
they describe a model chosen in advance.

Ranking itself is unaffected — every candidate is scored on the same samples
with the same degrees of freedom. If you need a number to quote for the winner,
select on one subset of samples and report on another.

### 2) Case/control traits

Phenotype type provided via `--pheno-name` is detected automatically. 
If it is binary trait add `--prevalence` to calcaulte liability-scale R².
ROC AUC is caclulated automatically for case/control traits.

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

### 3) Correlation matrix

It is also possible to produce correlation matrix for used PGSs.

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 3-5 \
         --give-corr-matrix --out result
```

| ID | PGS_A | PGS_B | PGS_C |
|---|---|---|---|
| PGS_A | 1 | -0.0248890537084 | -0.0589350562717 |
| PGS_B | -0.0248890537084 | 1 | -0.0752174055941 |
| PGS_C | -0.0589350562717 | -0.0752174055941 | 1 |

### 4) Evaluate only certain pairs

You can specify which phenotypes from `--pheno` should be used to evaluate each score. 
Define the corresponding score–phenotype pairs in the file passed to `--pgs-pheno-pairs`:

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

### 5) Evaluate single PGS against all of the phenotypes

Sometimes you have only single PGS and you need to evalaute it agains all of the phenotypes in `--pheno`.

```sh
adjutant --bfile demo --score weights.txt 1 2 4 header-read \
         --pheno pheno.txt --covar covar.txt --give-best-pgs --out result
```

| PGS_ID | PHENO | TYPE | N | N_CASE | R2_INCR | R2_INCR_L95 | P |
|---|---|---|---|---|---|---|---|
| PGS_B | HEIGHT | QT | 300 | NA | 0.171924 | 0.099872 | 7.77127e-14 |
| PGS_B | DISEASE | CC | 300 | 148 | 0.128126 | 0.0649802 | 7.7551e-09 |

### 6) Use a score matrix you already have

Sometimes you already have a `.sscore` from a previous run, or from plink2 or adjutant, and want to skip the scoring step.
`--sscore` takes any `.sscore` from plink2 or adjutant, and needs no genotypes:

```sh
adjutant --sscore result.sscore --pheno pheno.txt --covar covar.txt \
         --give-best-pgs --out ranked
```

### 7) Drop the columns your run cannot fill

`.best_pgs` always has 20 columns, eight of which describe binary traits. A run
where every phenotype is quantitative fills all eight with `NA`, so roughly two
fifths of the file carries nothing. `no-na` omits them:

```sh
adjutant --sscore result.sscore --pheno metabolites.txt --covar covar.txt \
         --give-best-pgs no-na --out ranked
```

| dropped | when |
|---|---|
| `N_CASE`, `OR`, `AUC_NULL`, `AUC_FULL`, `AUC_INCR` | no phenotype is case/control |
| `R2L_INCR`, `R2L_INCR_L95`, `R2L_INCR_U95` | that, or `--prevalence` is absent |

The rule is the configuration, not the data. A column is dropped because the run
cannot populate it, not because its values happened to come out `NA`, so the
column set follows from the command line alone and two runs over the same cohort
always produce the same layout. A column left empty because every fit failed to
converge is therefore kept.

One binary phenotype anywhere keeps the case/control columns for every row,
including the quantitative ones, where they stay `NA`. A column set that varied
from row to row would not be a table.

### 8) Pin the columns, and their positions

`no-na` decides the layout for you, which is convenient to read and awkward to
parse: the column set depends on the run. `cols=` inverts that — you name the
columns and they come out in the order you named them:

```sh
adjutant --sscore result.sscore --pheno pheno.txt --covar covar.txt \
         --give-best-pgs cols=PGS_ID,PHENO,R2_INCR,P --out ranked
```

| PGS_ID | PHENO | R2_INCR | P |
|---|---|---|---|
| PGS_B | HEIGHT | 0.171924 | 7.77127e-14 |

**This is how to make a positional parser safe.** `$3` is `R2_INCR` because you
put it there, and no future version of adjutant can move it — adding a column
changes the default layout, never a named one. Without `cols=`, index by name
off the header line instead; `.best_pgs` has lost and gained columns before and
will again.

An unknown or repeated name is a parse error, raised before any genotype is
read. The names are the header's, without the leading `#`:

```
PGS_ID PHENO TYPE N N_CASE BETA OR SE P R2_NULL R2_FULL
R2_INCR R2_INCR_L95 R2_INCR_U95 R2L_INCR R2L_INCR_L95 R2L_INCR_U95
AUC_NULL AUC_FULL AUC_INCR
```

`cols=` and `no-na` are alternatives; give one or the other.

---

## III. Harmonisation

### 1) Align a published PGS to your target .bim

A published PGS rarely uses your dataset's allele coding. `--pgs-model` aligns
it, reimplementing [`bigsnpr::snp_match()`](https://github.com/privefl/bigsnpr/blob/master/R/snp_match.R).

```sh
adjutant --bim target.bim --pgs-model published_pgs.txt 1 2 3 4 5 --out harmonised
```

The five numbers are the columns holding (1) chromosome, (2) position, (3) effect allele, (4) other allele and (5) weight.

```
--pgs-model: 40 source variants, 40 found in the .bim by position.
--pgs-model: 14 ambiguous (A/T, C/G) dropped; use --allow-ambiguous to keep them.
--pgs-model: 26 matched (8 strand-flipped, 10 reversed), 0 unmatched on alleles.
```

`harmonised.matched` is the output, with the target's coordinates and alleles, and the weight reoriented to the target's effect allele:
| CHR | POS | ID | A1 | A2 | BETA |
|---|---|---|---|---|---|
| 1 | 2 | snp2 | A | G | 0.4 |
| 1 | 3 | snp3 | T | C | 0.4 |
| 1 | 4 | snp4 | C | T | -0.1 |

It feeds straight back into `--score`, with `ID` in column 3, `A1` in 4 and `BETA` in 6:

```sh
adjutant --bfile demo --score harmonised.matched 3 4 6 header-read --out result
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

### 2) Matching on rsID instead of position

If the weight file and your target come from **different genome builds**,
`chr:pos` matching finds nothing. A leading `id` joins on the variant ID
instead. Sometimes it is useful if you use [PRS-CS](https://github.com/getian107/PRScs):

```sh
adjutant --bim target.bim --pgs-model pgs_rsid.txt id 1 4 5 6 --out harmonised
```

Note `id` flag, then the column numbers: 1 for ID, 4 for effect allele, 5 for other allele, 6 for weight.
Chromosome and position go unused.

> [!IMPORTANT]
> **This trades a build problem for an rsID-versioning one.** 
> Because rsIDs can be merged or retired, source and target files based on different dbSNP releases may fail to match (or may incorrectly match) variants. 
> Unlike a genome-build mismatch, this issue cannot be detected from genomic coordinates because coordinates are not used in the comparison.
> Adjutant therefore reports the match rate for every run.

### 3) Rescuing strand-ambiguous variants

A/T and C/G variants look the same on either strand. Given frequencies from a
**comparable ancestry** on both sides, `--resolve-ambiguous` resolves them.
GWAS summary statistics often carry a frequency column, and the target's allele frequencies can be read from a `.afreq` file. 
The two frequencies must be **both below 0.4 or both above 0.6** to reliably call the strand; otherwise the variant is dropped. 
The target's frequency is in column 6 of the weight file, and the source frequency is in the `.afreq` file.

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
refused rather than guessed.

### 4) Merging several models

Sometimes you have several published models for the same trait or just numerous of them, and you want to score them in one batch.

`--pgs-models-list` harmonises several published models against the **same**
target and merges them into one matrix - one column per model, one row per
variant matched by at least one of them.

Give it a text file listing the weight files, one per line. **Absolute paths
are the safe choice** — a relative path resolves against the directory you run
from, not the directory the list lives in:

```
/data/pgs/pgs_alpha.txt
/data/pgs/pgs_beta.txt
```

Relative paths work too, and `example/models.txt` uses them so it runs from
`example/`.

Blank lines and `#` comments are skipped. Each column takes its name from the
file's basename up to the first dot, so `PGS000123.txt.gz` becomes `PGS000123`;
two paths yielding the same name is an error.

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

A model that does not carry a variant gets `0`, not a missing value. 

**Adjutant harmonises models before merging.** Two models can name opposite
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
> exactly `4.0` (= 2 × 0.4 × 5), `r = 1.0` across 200 samples.
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

20 columns, the same set for both trait types.

| Column | Meaning |
|---|---|
| `PGS_ID`, `PHENO` | which model, which phenotype |
| `TYPE` | `QT` or `CC` |
| `N`, `N_CASE` | complete cases; cases among them (`CC` only) |
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
