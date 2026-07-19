# MethylForge

**A reproducible Nextflow + Python pipeline for epigenome (DNA-methylation) analysis — from short- or long-read methylation data to imprinting-disorder calls, episignature classification, differential-methylation regions, and branded clinical / research reports.**

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2%20%E2%89%A523.10-23aa62.svg)
![Python](https://img.shields.io/badge/CLI-Python-3776ab.svg)

MethylForge is the **epigenome** counterpart to [CallForge](https://github.com/aduton1000/callforge) (sequence variants) and [FamilyForge](https://github.com/aduton1000/familyforge) (interpretation). It is **organism- and disease-agnostic**: change inputs and config, never code.

## What it does

- **Any methylation input** — short-read **bisulfite / EM-seq** (Bismark or bwa-meth), **long-read native 5mC** (Nanopore `dorado`→`modkit`, PacBio HiFi), pre-called `bedMethyl`/`bedGraph`/Bismark `.cov`, or (optional) Illumina EPIC arrays. **Genome-wide or targeted/capture.**
- **Imprinting-disorder analysis** — methylation at a curated panel of imprinting control regions with loss/gain calls mapped to disorders (Prader-Willi, Angelman, Beckwith-Wiedemann, Silver-Russell, Temple/Kagami-Ogata, TNDM, pseudohypoparathyroidism), plus multi-locus imprinting disturbance (MLID).
- **Episignature classification** — reference-based genome-wide methylation classifier (nearest-centroid by default; scikit-learn SVM if installed) with a calibrated confidence and an explicit *inconclusive* zone.
- **Differential methylation** — case-vs-control DMP/DMR calling (Mann-Whitney + BH-FDR, region merging) with gene annotation.
- **Methylation QC** — bisulfite/EM conversion efficiency, coverage, global methylation by context, inferred sex, pass/fail gates.
- **Report engine** — a **GBCF/Noguchi-branded clinical report** (with the ACMG/AMP severity colour convention: Pathogenic=red, Likely-pathogenic=orange→red, VUS=amber, Benign=green) and a **modern research report** with full DMR detail — each self-contained HTML with embedded figures.
- **FamilyForge interoperability** — every run emits a `candidates.json`-compatible findings file and a tolerant findings TSV that FamilyForge can ingest for downstream interpretation/reporting.

## Install

**This Mac (Apple Silicon):**
```bash
git clone https://github.com/aduton1000/methylforge.git
cd methylforge
mamba env create -f env/methylforge.yml     # core analysis + report engine
conda activate methylforge
pip install -e .
methylforge version
```

**Lab Ubuntu HPC:** see [`docs/hpc_deployment.md`](docs/hpc_deployment.md) (conda + Singularity + a `hpc_slurm` profile and an `lmod` module so a run is one command for the whole lab).

> The heavy aligners (Bismark, bwa-meth, MethylDackel, minimap2/modkit) live in their own per-stage conda envs (`env/align_short.yml`, `env/align_long.yml`, `env/extract.yml`) and build on demand. On Apple Silicon create those under Rosetta (`CONDA_SUBDIR=osx-64 mamba env create -f env/align_short.yml`); Linux/HPC is native. ONT `dorado` basecalling (GPU) is done upstream and is not a conda package.

## Quickstart — prove it works on synthetic data

```bash
methylforge testdata -o mf_demo          # synthetic methylation for 3 samples
bash test/run_demo.sh mf_demo            # runs PWS / Sotos / normal end-to-end
# or the programmatic assertion suite:
python tests/selftest.py
```

This produces branded clinical + research HTML reports and FamilyForge exports, and demonstrates the expected calls: **PWS case → Prader-Willi (SNRPN hypermethylation)**, **Sotos case → positive episignature + a DMR**, **control → all normal**.

## Run on real data (Nextflow)

```bash
# short-read bisulfite/EM-seq or long-read, driven by a sample sheet
nextflow run main.nf -profile local,conda \
  --samplesheet samples.csv --reference GRCh38.fa \
  --episign_ref ref/episig_reference.tsv --episign_labels ref/episig_labels.tsv

# on the HPC
nextflow run main.nf -profile hpc_slurm,singularity --samplesheet samples.csv --reference GRCh38.fa
```

`samples.csv`: `sample,fastq_1,fastq_2,read_type,modbam,methylation,group` — `read_type` ∈ `short|long|precalled`.

See the full [**MANUAL**](MANUAL.md) for every command, stage, resource, and best-practice note.

## Pipeline

```
FASTQ / modBAM / bedMethyl / .cov / EPIC IDAT
        │
   ┌────┴─────────────────────────────┐
short-read                        long-read
  QC/trim (fastp/Trim Galore)      align (minimap2 -y)
  align (Bismark | bwa-meth)       pileup (modkit)
  dedup + extract (MethylDackel)        │
        └───────────► per-CpG methylation ◄── pre-called input
                              │
        ┌─────────────────────┼───────────────────────┐
     QC (stage4)     imprinting (5)   episignature (7)   DMR (6)
        └─────────────────────┼───────────────────────┘
                       report engine (9): clinical (GBCF-branded) + research
                              │
                    FamilyForge export (candidates.json + findings TSV)
```

## Important

MethylForge is research software. Panel coordinates, thresholds, and episignature classifiers **must be locally validated** — ideally against **ancestry-matched controls**, since methylation references are overwhelmingly European-ancestry — before any diagnostic use. Findings must be interpreted by a qualified professional with orthogonal confirmation where indicated.

## Citation / license

MIT licensed. If you use MethylForge, please cite this repository. © 2026 aduton1000.
