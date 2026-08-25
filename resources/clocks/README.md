# Epigenetic clocks for MethylForge

MethylForge's clock engine (`methylforge clock`, `run --clock`) is generic: it
consumes a per-clock coefficient TSV (`feature`, `coefficient`, optional `mean`)
plus a `.meta.json` transform spec. Coefficient **values are not bundled** — they
belong to their publications — so build them from an authoritative source with
`bin/prepare_clocks.R`, which writes MethylForge-format files into this directory.

## Build

```
Rscript bin/prepare_clocks.R            # -> _local/refs/clocks/ (gitignored data)
```
Coefficient *values* come from the `dnaMethyAge` package (authoritative, pure-R;
not fabricated) and are written to `_local/refs/clocks/` (reference data, kept out
of git per repo convention). Each clock is emitted twice: `<clock>.probe.tsv`
(EPIC/450K probe-keyed, for array betas or WGBS after `methylforge array`) and
`<clock>.hg38.tsv` (GRCh38 `chrom:pos`, for sequencing — the twin methylomes),
each with a `.meta.json` transform spec. Probe→hg38 mapping uses the 450K/EPIC
minfi annotation + UCSC liftOver (`bin/probe_hg38_map.R`), no sesame cache needed.
Options: `--out <dir>`, `--clocks id,id`, `--no-coords`, `--chain <chain>`.
This `manifest.tsv` (also written to the output dir) is the machine-readable
registry — clock id, transform, adult_age, CpG count, tissue, unit, dnaMethyAge
source id, and DOI.

## Transforms (`.meta.json`)

`transform: "horvath"` applies Horvath's calibrated age transform (below the
`adult_age` knot the predictor is anti-logged); `transform: "none"` reads the
linear predictor directly as years (or as a pace, for DunedinPACE).

## `manifest.tsv`

The registry of clocks this directory can build — clock id, transform, adult_age,
CpG count, tissue, unit, source id in `dnaMethyAge`, and DOI. Pediatric clocks are
flagged: they matter here because the study includes **infants**, where adult
blood clocks are extrapolating far outside their training range — prefer a
pediatric clock and always report cell-composition-adjusted (intrinsic)
acceleration (`clocks.intrinsic_acceleration`).

> Clocks are applied per sample, so they work at any n. But age acceleration in a
> tiny cohort is descriptive only, and is confounded by cell composition and by
> chronological-age differences — interpret alongside the deconvolution output.
