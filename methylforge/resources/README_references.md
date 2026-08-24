# Reference & coefficient file formats (v0.2 stages)

MethylForge does **not** redistribute cell-type references or clock coefficients
(they belong to their publications). Supply your own, validated for your platform
and genome build. All files are plain TSV and are keyed by **feature id**, which
must match your data's keyspace: `chrom:pos` for sequencing (WGBS/EM-seq), or a
probe id (`cg########`) for array betas. Map array-probe-keyed references to
coordinates for sequencing data with `methylforge array`.

## Cell-type deconvolution reference (`--reference`)
Feature x cell-type mean methylation (beta 0..1); first column = feature id.

```
feature        CD4T    CD8T    NK      Bcell   Mono    Neu
chr1:100234    0.82    0.80    0.10    0.75    0.12    0.11
chr1:100512    0.05    0.06    0.55    0.04    0.61    0.60
...
```
Recommended sources: Salas et al. 2022 12-cell IDOL panel (blood); IDOL 450K/EPIC
libraries. For WGBS/EM-seq, lift the panel's probes to `chrom:pos` on your build.
Estimators: `CP` (constrained projection, default) or `RPC` (robust).

## Epigenetic-clock coefficients (`--clock name=path`)
Columns: `feature`, `coefficient`, optional `mean` (for imputing markers missing
from the sample). A special `(Intercept)` row carries the intercept.

```
feature        coefficient   mean
(Intercept)    0.6955        NA
cg00075967     0.1290        0.53
chr6:31234567  -0.0844       0.61
...
```
Non-linear transforms (e.g. the Horvath age transform) go in a sidecar
`<clock>.meta.json` (or pass `meta_path`):

```json
{"name": "Horvath2013", "transform": "horvath", "adult_age": 20}
```
`transform` is `none` (additive clocks: Hannum, PhenoAge, DunedinPACE-style) or
`horvath`. Clocks are per-sample; use `--chrono-age` for delta age, and pass
deconvolution proportions to `clocks.intrinsic_acceleration()` for the
composition-independent (intrinsic) signal.

## Allele-specific methylation input (`asm`)
Either two haplotype methylation files (`--hp1/--hp2`, any supported per-cytosine
format — e.g. `modkit pileup --partition-tag HP`), or a long allele table
(`--allele-table`) with columns `chrom,pos,allele,meth,cov`. `--candidate <file>`
runs the single-sample heuristic scan (flags loci to confirm by phasing).
