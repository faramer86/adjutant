# Example data

The dataset every command in the top-level [README](../README.md) runs against.
300 samples, 40 variants — small enough that real output fits in a table, big
enough for the statistics to mean something.

| File | What it is |
|---|---|
| `demo.bed` / `.bim` / `.fam` | genotypes, PLINK 1 binary form (used by the examples) |
| `demo.pgen` / `.pvar` / `.psam` | the same genotypes, PLINK 2 form |
| `demo.afreq` | allele frequencies, for `--read-freq` |
| `weights.txt` | three candidate models, `PGS_A` / `PGS_B` / `PGS_C` |
| `pheno.txt` | `HEIGHT` (quantitative) and `DISEASE` (case/control) |
| `covar.txt` | `AGE` and `PC1`, unrelated to the signal |
| `partitions.txt` | two variant subsets, for `--score-partitions` |
| `pairs.txt` | phenotype/PGS pairs, for `--pgs-pheno-pairs` |
| `target.bim` | a harmonisation target |
| `published_pgs.txt` | a "published" PGS whose alleles are deliberately mis-oriented |
| `pgs_alpha.txt` / `pgs_beta.txt` | two overlapping models, for `--pgs-models-list` |
| `models.txt` | the list file naming those two, one path per line |
| `pgs_rsid.txt` | an rsID-keyed PGS whose positions are a different build |
| `mhc.bed` / `.bim` / `.fam` | 24 variants straddling the MHC boundaries, for `--exclude-mhc-*` |
| `mhc_weights.txt` | weights over all 24, for `--score` |
| `mhc_pgs.txt` | the same variants keyed on chr/pos, for `--pgs-model` |

## The MHC fileset is positioned, not random

`demo` is entirely chromosome 1 at positions 0–39, so the `--exclude-mhc-*`
flags remove nothing from it. `mhc.*` exists to exercise them: 200 samples and
24 variants placed *exactly* on the interval edges — one at each `start` and
`end`, one at each `start-1` and `end+1`, for all four regions — plus a
chromosome 7 variant at an MHC position, so a filter that forgot to check the
chromosome fails visibly.

| flag | excludes |
|---|---|
| `--exclude-mhc-hg19` | 10 |
| `--exclude-mhc-hg38` | 10 |
| `--exclude-xmhc-hg19` | 9 |
| `--exclude-xmhc-hg38` | 13 |

The two `mhc` builds excluding the same **count** while disagreeing about four
variants is deliberate: it means a test comparing only counts would pass with
the hg19 and hg38 coordinates swapped.

Four of the alleles are A/T or C/G, so a `--pgs-model` run against this fileset
also exercises the ambiguous-variant path rather than matching cleanly.

## It is not random noise

`HEIGHT` is **generated from `PGS_B`**, so `--give-best-pgs` has a real answer
to find: it ranks `PGS_B` first at R² ≈ 0.17 while the two decoys get an
`R2_INCR_L95` of exactly 0, and under `--folds 5` the decoys go negative.
`DISEASE` splits on the same signal, tuned to a realistic OR ≈ 2.4 / AUC ≈ 0.70
rather than a near-perfect separation that would teach nothing.

`published_pgs.txt` has its alleles deliberately swapped and strand-flipped, so
harmonisation reports a mix of `EXACT`, `REVERSED`, `FLIPPED` and `AMBIGUOUS`
rather than a clean pass.

`pgs_alpha.txt` and `pgs_beta.txt` cover rows 1–25 and 15–40 of `target.bim`, so
`--pgs-models-list` has to zero-fill in both directions. Every third `beta` row
names the **opposite** effect allele, which is what makes the merge report
`REVERSED` rows instead of taking the weights at face value. Their 11 overlapping
source rows collapse to only 3 nonzero in both columns, because the rest are
strand-ambiguous — a worked example of why a `0` in `.matched` needs
`.match_report` to interpret.

`pgs_rsid.txt` is keyed on rsID and its positions are deliberately shifted by
1,000,000 — a stand-in for a weight file published against a different genome
build. Matched on `chr:pos` it finds **0 of 40** variants; matched on ID it finds
**40 of 40**. That contrast is the argument for `--pgs-model ... id`, shown
rather than asserted.

## Provenance

The genotypes are simulated, not real people:

```sh
plink2 --dummy 300 40 acgt --make-pgen --seed 20260803 --out demo
plink2 --pfile demo --make-bed --out demo
plink2 --pfile demo --freq   --out demo

plink2 --dummy 200 24 acgt --make-bed --seed 20260805 --out mhc
# mhc.bim is then rewritten with the chosen chromosomes, positions and IDs;
# the .bed holds genotype codes indexed by row, so relabelling rows is valid.
```

Everything else is derived from those genotypes deterministically — the "noise"
in the phenotypes is a function of the sample index, not `rand()` — which is
what lets the README quote exact numbers.
