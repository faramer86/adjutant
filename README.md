# Adjutant 🫡

**PLINK 2, cut down to polygenic scores** — computing them, evaluating them, and
aligning published ones to your data.

One binary, no runtime dependencies, and `--score` output that is byte-identical
to stock plink2.

> **Most of this code is [PLINK 2.0](https://github.com/chrchang/plink-ng)'s**,
> and the credit belongs there. `adjutant` was made by deleting commands and reusing functions, not by
> writing a scorer. If you publish with it, hit a star for this repository and cite PLINK 2 — see [NOTICE](NOTICE).

---

## Three modes

| | Input | Output | Flags |
|---|---|---|---|
| **Score** | genotypes + weights | `.sscore` | `--score` |
| **Evaluate** | scores + a phenotype | `.best_pgs`, `.pgs_corr` | `--give-best-pgs`, `--give-corr-matrix` |
| **Harmonise** | weights + a `.bim` | `.matched` | `--pgs-model` |

Score and evaluate run in one command. Evaluate also takes a score matrix you
already have (`--sscore`). Harmonise is standalone — two text files in, two out,
no genotypes.

---

## Build

```sh
make release          # -> bin/adjutant
```

Needs a C++17 compiler, zlib, and a BLAS/LAPACK (Accelerate on macOS, OpenBLAS
elsewhere). Everything else is vendored.

Pre-built macOS arm64 and Linux x86-64 binaries — including a fully static one
that needs no glibc — are attached to the GitHub release rather than committed
here.

---

## Quick start

**Every example below runs against [`example/`](example/)** — a complete,
committed set of working input files. Clone the repo and they run as written:

```sh
cd example
```

300 samples and 40 variants, small enough to read, and not random noise:
`HEIGHT` is genuinely generated from `PGS_B`, so the tool has a real answer to
find. The genotypes ship as both a `.bed/.bim/.fam` trio and a
`.pgen/.pvar/.psam` one; the examples use `--bfile`, and `--pfile demo` works
identically.

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

Copy any of them as a template for your own data — see
[example/README.md](example/README.md) for how they were made.

---

## Scoring

### One model

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --out result
```

`1 2` are the columns holding the variant ID and the effect allele;
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

The `FID` column comes from the `.fam`, which always carries one — `--pfile`
output on a `.psam` without FIDs has no such column. Everything else is
identical between the two.

### Split a model into pathways

`--score-partitions` takes a 0/1 mask per variant and scores each subset
separately, in the same pass:

```
SNP    PATHWAY_A  PATHWAY_B
snp0   1          1
snp1   0          1
```

```sh
adjutant --bfile demo --score weights.txt 1 2 header-read --score-col-nums 4 \
         --score-partitions partitions.txt --out result
```

| FID | IID | ALLELE_CT | NAMED_ALLELE_DOSAGE_SUM | PGS_B_AVG | PGS_B_PATHWAY_A_AVG | PGS_B_PATHWAY_B_AVG |
|---|---|---|---|---|---|---|
| 0 | per0 | 80 | 34 | 0.01375 | 0.0075 | 0.005 |
| 0 | per1 | 80 | 34 | 0.0375 | 0.04 | 0.0525 |
| 0 | per2 | 80 | 48 | 0.02125 | -0.0025 | 0.045 |

`PGS_B_AVG` is the unpartitioned total. Each partition divides by **its own**
variants' allele count, so partitions of different sizes stay comparable — and
therefore do **not** sum to the total. Use `cols=+scoresums` if you need sums,
which do.

### Drop the MHC

The MHC's long-range LD and extreme polymorphism make it behave unlike the rest
of the genome, so PGS work often excludes it. Pick the flag matching **your
data's build** — there is no auto-detection, and the wrong one silently drops
the wrong window:

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

Every run logs the interval it applied, so a log says which coordinates were
used rather than leaving you to infer them from the flag name.

The flags work the same way in harmonisation, where they apply to **both** the
`.bim` and the weight file, so excluded variants are absent from the
`.match_report` rather than reported as unmatched:

```sh
adjutant --bim mhc.bim --pgs-model mhc_pgs.txt 1 2 3 4 5 --exclude-xmhc-hg38 --out h
```

> [!IMPORTANT]
> **The extended region does not contain the plain one.** The two come from
> different published definitions: the GRC MHC region ends **79,933 bp past** the
> xMHC end, in both builds. So `--exclude-xmhc-hg19` reaches 2.75 Mb further 5′
> but stops short 3′, and can exclude **fewer** variants than
> `--exclude-mhc-hg19` — on the `example/mhc` fileset it excludes 9 against 10.
>
> This is intrinsic to the source definitions, not a quirk of this tool; no
> pairing of published MHC and xMHC coordinates nests. adjutant reports each
> region exactly as published rather than silently widening one to fit the other.
> If you want strictly everything, run the plain flag as well in a separate pass
> and intersect the results.

---

## Evaluation

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
interval before you read the ranking.

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

`--folds 1` is the default and is bit-identical to omitting the flag.

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

## Harmonisation

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

`harmonised.matched` feeds straight back into `--score`:

| ID | A1 | BETA |
|---|---|---|
| snp2 | A | 0.4 |
| snp3 | T | 0.4 |
| snp4 | C | -0.1 |

`harmonised.match_report` keeps every variant, including the dropped ones:

| STATUS | count | meaning |
|---|---|---|
| `EXACT` | 8 | alleles already match |
| `REVERSED` | 10 | effect/other swapped — **beta negated** |
| `FLIPPED` | 8 | opposite strand |
| `AMBIGUOUS` | 14 | A/T or C/G — strand undecidable, dropped |

**Weights must be betas.** Negating an odds ratio is meaningless; adjutant warns
if every weight is positive with mean in [0.5, 2.0].

### Rescuing strand-ambiguous variants

A/T and C/G variants look the same on either strand. Given frequencies from a
**comparable ancestry** on both sides, `--resolve-ambiguous` decides them:

```sh
adjutant --bim target.bim --pgs-model published_pgs.txt 1 2 3 4 5 6 \
         --read-freq demo.afreq --resolve-ambiguous --out resolved
```

```
--pgs-model: 14 ambiguous; 11 resolved by frequency (MAF < 0.4 on both sides),
--pgs-model: 37 matched (8 strand-flipped, 13 reversed), 0 unmatched on alleles.
```

26 matched becomes 37. The 3 left over have a MAF too near 0.5 to call and are
refused rather than guessed. Column 6 is the source effect-allele frequency.

### Merging several models

`--pgs-models-list` takes a file of weight-file paths and merges them into one
matrix, zero-filling variants a model does not carry. Each column is named from
that file's basename.

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
contract, and it is why the build uses `-ffp-contract=off` — the flag stops the
compiler fusing multiply-add in the scoring dot product, which is what keeps
`-O3` identical to `-O2` and to upstream.

You do not have to take that on trust. If you have plink2 installed, check it
against your own data:

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
