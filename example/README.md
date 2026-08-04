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

## It is not random noise

`HEIGHT` is **generated from `PGS_B`**, so `--give-best-pgs` has a real answer
to find: it ranks `PGS_B` first at R² ≈ 0.17 while the two decoys get an
`R2_INCR_L95` of exactly 0, and under `--folds 5` the decoys go negative.
`DISEASE` splits on the same signal, tuned to a realistic OR ≈ 2.4 / AUC ≈ 0.70
rather than a near-perfect separation that would teach nothing.

`published_pgs.txt` has its alleles deliberately swapped and strand-flipped, so
harmonisation reports a mix of `EXACT`, `REVERSED`, `FLIPPED` and `AMBIGUOUS`
rather than a clean pass.

## Provenance

The genotypes are simulated, not real people:

```sh
plink2 --dummy 300 40 acgt --make-pgen --seed 20260803 --out demo
plink2 --pfile demo --make-bed --out demo
plink2 --pfile demo --freq   --out demo
```

Everything else is derived from those genotypes deterministically — the "noise"
in the phenotypes is a function of the sample index, not `rand()` — which is
what lets the README quote exact numbers.
