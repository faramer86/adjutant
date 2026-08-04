# Test data

Everything the four verification suites read. All of it is committed so that
`./score_check.sh`, `./input_matrix.sh`, `./pgs_eval_check.sh` and
`./demo_give_best_pgs.sh` run on a fresh clone with nothing else to download.

| Path | Used by | What it is |
|---|---|---|
| `TEST_PHASED_VCF/` | `score_check.sh`, `pgs_eval_check.sh`, `demo_give_best_pgs.sh` | 2504 samples, chrM and chr21 subsets |
| `input_matrix/` | `input_matrix.sh` | the same cohort as both a `.pgen` and a `.bed` fileset, plus weight files |
| `pgs_eval_fixtures/` | `pgs_eval_check.sh` | reference values from statsmodels and R |
| `pgs_match_fixtures/` | `pgs_eval_check.sh` | reference output from `bigsnpr::snp_match()` |
| `score_fixtures/` | `score_check.sh` | weight files for the golden cases |

## Provenance of the genotypes

The genotype data is from the **1000 Genomes Project** (phase 3, and an hg38
lift-over), identifiable by its `HG#####` / `NA#####` sample IDs.

The 1000 Genomes panel is open-access, published under
[full consent for public release](https://www.internationalgenome.org/),
including redistribution — which is why these files can sit in a public
repository. It carries no restricted or identifiable clinical information.

Phenotypes in the suites are **simulated**, not real: they are generated from
the genotypes by the test scripts themselves, so nothing here associates a real
trait with a real person.

## Reference fixtures

`pgs_eval_fixtures/` and `pgs_match_fixtures/` hold values computed once by
statsmodels, R and `bigsnpr`, then committed. Neither Python nor R is needed to
run the suites — that is the point of committing them. The generator scripts sit
alongside the values they produce.

## Golden `--score` output

`../score_regression_golden/` holds 126 files captured from stock plink2 and is
the contract `score_check.sh` enforces: adjutant's `--score` output must match
them byte for byte. A run compares 89 of them — `gone_*` (proof that a removed
command is refused) and `_failures.txt` (the expected-nonzero-exit list) are
excluded from the diff by design.

The `.log` files there are **normalised** by `score_regression.sh` before being
compared or committed — version, hostname, working directory, seed and RAM are
replaced with `<VERSION>`, `<HOST>`, `<CWD>`, `<SEED>`, `<MEM>`. Those fields
differ on every machine, so a raw log would both fail the comparison and leak
the path it was generated in.
