# Test data

Reference data and fixtures for verifying a build. `make pgs_eval_test` reads
the fixture directories directly and needs nothing else; the golden `--score`
output under `../score_regression_golden/` is there to diff a build against by
hand.

| Path | What it is |
|---|---|
| `TEST_PHASED_VCF/` | 2504 samples, chrM and chr21 subsets |
| `input_matrix/` | the same cohort as both a `.pgen` and a `.bed` fileset, plus weight files |
| `pgs_eval_fixtures/` | reference values from statsmodels and R, read by `pgs_eval_test` |
| `pgs_match_fixtures/` | reference output from `bigsnpr::snp_match()` |
| `score_fixtures/` | weight files for the golden `--score` cases |

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
the contract adjutant's `--score` output must match byte for byte. 89 are
directly comparable; the rest are `gone_*` (captured from commands this build
removes, kept as a record of what they used to produce) and `_failures.txt`,
the list of cases expected to exit nonzero.

The `.log` files there are **normalised** — version, hostname, working
directory, seed and RAM are replaced with `<VERSION>`, `<HOST>`, `<CWD>`,
`<SEED>`, `<MEM>`. Those fields differ on every machine, so a raw log would
both defeat the comparison and leak the path it was generated in. Each log's
"Options in effect" block records the exact command that produced its
`.sscore`, which is what makes the goldens reproducible by hand.
