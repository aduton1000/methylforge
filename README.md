# MethylForge

**A reproducible Nextflow + Python pipeline for epigenome (DNA-methylation) analysis — from short- or long-read methylation data (or arrays) to imprinting-disorder calls, episignature classification, differential-methylation regions, and branded clinical / research reports.**

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2%20%E2%89%A523.10-23aa62.svg)
![Python](https://img.shields.io/badge/CLI-Python%20%E2%89%A53.10-3776ab.svg)

MethylForge is the **epigenome** counterpart to [CallForge](https://github.com/aduton1000/callforge) (sequence variants) and [FamilyForge](https://github.com/aduton1000/familyforge) (interpretation). It is **organism- and disease-agnostic**: change inputs and config, never code. Every input — bisulfite/EM-seq, long-read native 5mC, pre-called methylation, or Illumina EPIC/450K arrays — is normalised to a single per-cytosine schema, so all downstream analysis is caller-agnostic.

## The biology it addresses

DNA methylation (a methyl group on cytosine, almost always at **CpG** sites) is a stable epigenetic mark controlling imprinting, X-inactivation and gene expression. MethylForge quantifies it as a per-cytosine **beta value** (methylated fraction, 0–1) and answers three clinically important questions: is methylation aberrant at an **imprinting control region** (→ Prader-Willi, Angelman, Beckwith-Wiedemann, Silver-Russell, Temple/Kagami-Ogata, transient neonatal diabetes, pseudohypoparathyroidism)? does the genome-wide profile match a disorder **episignature**? and where does methylation **differ** between groups (DMR)?

## What it does

- **Any methylation input** — short-read **bisulfite / EM-seq** (Bismark or bwa-meth), **long-read native 5mC** (Nanopore `dorado`→`modkit`, PacBio HiFi), pre-called `bedMethyl`/`bedGraph`/Bismark `.cov`, or **Illumina EPIC/450K** arrays (minfi). Genome-wide or targeted.
- **Methylation QC** — bisulfite/EM conversion efficiency (from lambda/pUC19 spike-ins, a provided value, or CHH), coverage, global methylation by context, inferred sex, pass/fail gates.
- **Imprinting-disorder analysis** — methylation at a curated panel of imprinting control regions with loss/gain calls mapped to disorders, plus a multi-locus imprinting disturbance (MLID) flag.
- **Episignature classification** — reference-based genome-wide classifier (nearest-centroid by default; scikit-learn SVM if installed) with calibrated confidence and an explicit *inconclusive* zone; includes a reference **builder + leave-one-out validator**.
- **Differential methylation** — case-vs-control DMP/DMR calling (Mann-Whitney / Welch-t on M-values + BH-FDR, region merging) with gene annotation.
- **Report engine** — a GBCF/Noguchi-branded **clinical** report (ACMG/AMP severity colours) and a modern **research** report (full DMR detail), each self-contained HTML, plus `results.json` and an **Excel findings workbook**.
- **FamilyForge interoperability** — every run emits a `candidates.json`-compatible findings file and a tolerant TSV that FamilyForge can ingest for downstream interpretation/reporting.

## Pipeline architecture

Two layers: a **Nextflow DSL2** layer (per-sample alignment/extraction + orchestration) and a **pure-Python** analysis + report engine (runs anywhere). Three modes selected with `--mode`: `sequencing` (default), `array`, `cohort_qc`.

| # | Stage (Nextflow module) | Tool(s) | Standalone CLI |
|:--|:--|:--|:--|
| 0 | Array ingest (`ARRAY_INGEST`/`ARRAY_TO_METHYL`) | minfi/sesame, `epic_process.R` | `methylforge array` |
| 1 | Read QC & trimming (`QC_TRIM`) | FastQC, fastp | — |
| 2 | Alignment — short (`ALIGN_SHORT`) / long (`ALIGN_LONG`) | Bismark \| bwa-meth · minimap2+modkit | — |
| 3 | Methylation extraction (`EXTRACT`) | MethylDackel | — |
| 4 | Methylation QC | `qc.py` | `methylforge qc` |
| 5 | Imprinted-DMR analysis | `imprinting.py` | `methylforge imprinting` |
| 6 | Differential methylation | `cohort.py`+`dmr.py` | `methylforge cohort`→`dmr` |
| 7 | Episignature classification | `episignature.py` | `methylforge episign` / `build-reference` |
| 8 | Region annotation | `annotate.py` | (via `--genes-bed`) |
| 9 | Reports + hand-off (`ANALYSE_REPORT`) | `report.py`, `export_excel.py`, `export_familyforge.py` | `methylforge run` / `findings-xlsx` |
| — | Batch cohort QC (`COHORT_QC`) | `dashboard.py` | `methylforge cohort-qc` |

The short-read (QC→align→extract), long-read (align→pileup), pre-called, and array paths all converge on the normalised per-CpG methylation, which feeds the shared analysis engine. See the [architecture diagram](docs/manual/figures/architecture.png).

## Install

**Local (Linux or macOS / Apple Silicon):**
```bash
git clone https://github.com/aduton1000/methylforge.git
cd methylforge
mamba env create -f env/methylforge.yml     # core analysis + report engine (pure Python)
conda activate methylforge
pip install -e .                            # run from the repo root
methylforge version
```
That covers everything operating on already-called methylation. To start from **raw reads / IDATs**, also create the per-stage envs (`env/align_short.yml`, `env/align_long.yml`, `env/extract.yml`, `env/array.yml`).

> On **Apple Silicon** the core env is native (arm64); the heavy-aligner envs are bioconda `osx-64` only — create them under Rosetta: `CONDA_SUBDIR=osx-64 mamba env create -f env/align_short.yml`. ONT `dorado` basecalling (GPU) is done upstream and is not a conda package.

**Lab HPC (SLURM + Singularity/Apptainer):** see [`docs/hpc_deployment.md`](docs/hpc_deployment.md).

## Quickstart — prove it works on synthetic data

```bash
methylforge testdata -o mf_demo          # synthetic methylation for 3 samples (+ cohort + array)
bash test/run_demo.sh mf_demo            # PWS / Sotos / control end-to-end
python tests/selftest.py                 # programmatic assertions
```
Expected: **PWS case → Prader-Willi (SNRPN hypermethylation)**, **Sotos case → positive episignature + a DMR**, **control → all normal**.

## Run on real data (Nextflow)

The classic-DSL2 modules use Nextflow's v1 language parser — `export NXF_SYNTAX_PARSER=v1` once per shell. Profiles are **composable**: one execution profile (`local` | `hpc_slurm`) + one packaging profile (`conda` | `mamba` | `docker` | `singularity` | `apptainer`).

```bash
export NXF_SYNTAX_PARSER=v1

# sequencing (default): short/long/precalled -> analysis + reports + auto batch QC
nextflow run main.nf -profile local,conda --mode sequencing \
  --samplesheet samples.csv --reference /path/to/GRCh38.fa \
  --episign_ref /path/to/episig_reference.tsv --episign_labels /path/to/episig_labels.tsv

# array: Illumina EPIC/450K IDAT -> beta (minfi) -> analysis + reports
nextflow run main.nf -profile local,conda --mode array \
  --idat_dir /path/to/idats --array_sheet /path/to/minfi_sheet.csv

# cohort_qc: batch QC dashboard only, from a samplesheet of methylation files
nextflow run main.nf -profile local,conda --mode cohort_qc --samplesheet /path/to/methyl_sheet.csv

# on the HPC
nextflow run main.nf -profile hpc_slurm,singularity --mode sequencing \
  --samplesheet samples.csv --reference /path/to/GRCh38.fa
```

**Sample sheet** (`sequencing` mode): `sample,fastq_1,fastq_2,read_type,modbam,methylation,group` — `read_type` ∈ `short|long|precalled` (see `test/samplesheet_precalled.csv`). Every `sequencing`/`array` run also emits a **batch QC dashboard** (`cohort_qc.html`); disable with `--cohort_qc false`.

## Outputs

Per sample under `--outdir`: `qc/`, `align/`, `methylation/`, and `report/` with the **clinical** and **research** HTML reports, `<sample>.results.json`, `<sample>.findings.xlsx`, and the FamilyForge hand-off (`<sample>.familyforge.json` + `epi_findings.tsv`). Per run: `cohort_qc.html` + `cohort_qc_summary.tsv` and `pipeline_info/` (Nextflow execution report/timeline/trace).

## Reference data

- **Imprinted-DMR panel** — bundled `methylforge/resources/imprinted_dmrs.tsv` (`--imprint_panel` to override); **validate coordinates/thresholds locally**.
- **Episignature reference** — build + validate with `methylforge build-reference` (`--episign_ref`/`--episign_labels`); use ancestry-matched controls.
- **Spike-ins** — append lambda + pUC19 to the reference with `bin/prepare_reference.sh` to measure conversion from the data.
- **EPIC manifest** — supplied by you (licensed); `bin/epic_process.R` derives coordinates via minfi/sesame (+ optional hg19→hg38 liftOver via `--array_chain`).

## Documentation

The full **[User Manual](docs/manual/methylforge_manual.pdf)** ([Markdown source](docs/manual/methylforge_manual.md)) documents every concept, stage, input schema, reference, configuration knob, and troubleshooting note. Build it with `pandoc metadata.yaml methylforge_manual.md -o methylforge_manual.pdf --pdf-engine=xelatex --toc --number-sections`.

## Important

MethylForge is **research software**. Panel coordinates, QC thresholds, and episignature classifiers **must be locally validated** — ideally against **ancestry-matched controls**, since public methylation references are overwhelmingly European-ancestry — before any diagnostic use. Findings must be interpreted by a qualified professional with orthogonal confirmation where indicated.

## Citation / license

MIT licensed (see [`LICENSE`](LICENSE)). If you use MethylForge, please cite it (see [`CITATION.cff`](CITATION.cff)). © 2026 Bright Adu, Genomics & Bioinformatics Core, Noguchi Memorial Institute for Medical Research, University of Ghana, Legon, Ghana.
