# MethylForge — User Manual

MethylForge is a reproducible epigenome (DNA-methylation) analysis pipeline. It takes methylation data from short- or long-read sequencing (genome-wide or targeted), performs QC, imprinting-disorder analysis, episignature classification, and differential-methylation calling, and produces branded clinical and research reports plus a FamilyForge-compatible export. This manual covers installation, inputs, every command and stage, the report engine, interoperability, and the standards/limitations you must observe.

---

## 1. Architecture

MethylForge has two layers:

- **A pure-Python analysis + report engine** (`methylforge/`) that depends only on `numpy/pandas/scipy/matplotlib/jinja2`. It runs anywhere (Mac, Linux, no heavy bioinformatics tools required) and is the part that produces QC, imprinting calls, episignatures, DMRs, and reports. This is the layer proven end-to-end by `tests/selftest.py`.
- **A Nextflow DSL2 wrapper** (`main.nf`, `modules/`, `subworkflows/`) that runs the heavy upstream stages (alignment, methylation extraction) with the standard tools and then calls the Python engine. It mirrors CallForge: two composable profile groups — execution (`local`/`hpc_slurm`) × packaging (`conda`/`mamba`/`docker`/`singularity`/`apptainer`).

```
methylforge/
  methylforge/            Python package (engine)
    io, qc, imprinting, dmr, episignature, annotate, plots, report, pipeline, cli, testdata
    export_familyforge    FamilyForge interoperability
    resources/            imprinted_dmrs.tsv, assets/ (GBCF + Noguchi logos)
    report_templates/     clinical_gbcf.html.j2 (branded), report.html.j2 (research)
  main.nf, nextflow.config, conf/, modules/, subworkflows/, env/
  test/, tests/, docs/
```

---

## 2. Installation

### 2.1 Mac (Apple Silicon)
```bash
git clone https://github.com/aduton1000/methylforge.git && cd methylforge
mamba env create -f env/methylforge.yml
conda activate methylforge
pip install -e .
methylforge version
python tests/selftest.py            # should print: self-test PASSED
```
The core env solves natively on `osx-arm64`. The **heavy-aligner** envs (`align_short.yml`, `align_long.yml`, `extract.yml`) contain tools without arm64 builds; create them under Rosetta:
```bash
CONDA_SUBDIR=osx-64 mamba env create -f env/align_short.yml
```

### 2.2 Linux HPC
See [`docs/hpc_deployment.md`](docs/hpc_deployment.md). Everything is native on linux-64; Singularity images are recommended for multi-user reproducibility.

### 2.3 Optional back-ends
- `pip install -e ".[ml]"` adds scikit-learn → the episignature classifier uses a linear SVM instead of the nearest-centroid fallback (automatic).
- `env/array.yml` (R/Bioconductor `minfi`/`sesame`) is only needed for Illumina EPIC-array input.

---

## 3. Inputs

### 3.1 Methylation data (any of)
| Input | Extension | Produced by |
|---|---|---|
| Bismark coverage | `.cov` / `.cov.gz` | Bismark methylation extractor |
| bedGraph | `.bedGraph` | MethylDackel |
| ENCODE bedMethyl | `.bedmethyl` | `modkit` (ONT), PacBio |
| MethylForge / generic TSV | `.tsv` | any; header with `chrom,start,beta` (or `meth,cov`) |

All are normalised internally to `chrom, start, end, context, meth, cov, beta`.

### 3.2 Sample sheet (Nextflow)
CSV with header `sample,fastq_1,fastq_2,read_type,modbam,methylation,group`:
- `read_type` = `short` (bisulfite/EM FASTQs), `long` (modBAM/FASTQ), or `precalled` (a methylation file).
- `group` = `case`/`control` (only used for cohort DMR).

### 3.3 Reference resources
- **Imprinted-DMR panel** (`methylforge/resources/imprinted_dmrs.tsv`): `locus, gene, chrom, start, end, normal_meth, hyper_disorder, hypo_disorder, notes`. Edit/extend freely; **validate coordinates and thresholds locally before diagnostic use.**
- **Episignature reference** (TSV): first column `feature` (CpG id `chrom:pos`), one column per reference sample (beta 0–1). Plus a labels TSV `sample,class`.
- **Genes BED** (optional): `chrom,start,end,name` for DMR annotation.

---

## 4. Command-line reference

### `methylforge run` — end-to-end single sample
```bash
methylforge run SAMPLE.cov -s SAMPLE --family FAM01 \
  --episign-ref ref/episig_reference.tsv --episign-labels ref/episig_labels.tsv \
  --dmr-matrix cohort/matrix.tsv --dmr-groups cohort/groups.tsv --genes-bed genes.bed \
  -o out/
```
Outputs in `out/`: `SAMPLE.clinical.report.html` (GBCF-branded), `SAMPLE.research.report.html`, `SAMPLE.results.json`, `SAMPLE.familyforge.json`, `SAMPLE.epi_findings.tsv`.

### `methylforge qc` — methylation QC only
```bash
methylforge qc SAMPLE.cov
```

### `methylforge imprinting` — imprinted-DMR analysis only
```bash
methylforge imprinting SAMPLE.cov [--panel my_panel.tsv]
```

### `methylforge cohort` — build a feature × sample matrix from many samples
Turns a set of per-sample methylation files into the `feature × sample` beta matrix
that the DMR and episignature stages consume.
```bash
# from a samplesheet (sample,methylation[,group]) — also writes groups.tsv
methylforge cohort --samplesheet cohort_samplesheet.csv --out matrix.tsv --groups-out groups.tsv
# or from globs; restrict to episignature features or a targeted region
methylforge cohort --inputs 'meth/*.cov' --sites episig_features.txt --out episig_matrix.tsv
methylforge cohort --samplesheet ss.csv --region-bed panel.bed --out targeted_matrix.tsv --groups-out groups.tsv
```
Options: `--min-cov` (per-site coverage filter, default 5), `--min-samples-frac` (keep a
site only if covered in this fraction of samples, default 0.75), `--sites` (feature list
`chrom:pos`), `--region-bed`, `--context` (default CpG). The output feeds straight into
`methylforge dmr` (with the groups TSV) or episignature reference building.

### `methylforge cohort-qc` — MultiQC-style batch QC dashboard
Aggregates per-sample QC across a whole run into one self-contained HTML: a batch summary,
overview plots (conversion, coverage, callable CpGs, global methylation, spike-in controls,
pass/fail), and a **sortable** per-sample table with pass/fail badges and red-flagged failing
cells.
```bash
# run QC on many samples from a samplesheet (or globs), and/or aggregate existing results.json
methylforge cohort-qc --samplesheet cohort_samplesheet.csv --out cohort_qc.html --summary-out cohort_qc.tsv
methylforge cohort-qc --inputs 'meth/*.cov' --out cohort_qc.html
methylforge cohort-qc --results 'out/*/*.results.json' --out cohort_qc.html
```
Samples that fail any gate (conversion, coverage, callable CpGs, pUC19 over-conversion) are
flagged **REVIEW** with the offending cell highlighted, so a whole batch is eyeballed at once.

### `methylforge dmr` — cohort differential methylation
```bash
methylforge dmr --matrix matrix.tsv --groups groups.tsv [--genes-bed genes.bed]
```
DMP test selection is automatic: Welch t-test on M-values for small cohorts (<12 per group,
where the Mann-Whitney p-floor is too high to survive FDR) and the Mann-Whitney rank test
for larger cohorts.

### `methylforge episign` — episignature classification
```bash
methylforge episign --reference ref.tsv --labels labels.tsv --query query.tsv
```

### `methylforge build-reference` — build + validate an episignature reference
Turns a labelled cohort matrix (feature × sample beta + `sample,class` labels, with a
`Control` class) into a curated, cross-validated reference you can ship and reuse.
```bash
methylforge build-reference --matrix cohort_matrix.tsv --labels labels.tsv \
  --features-per-class 40 --min-delta 0.15 -o my_reference/
```
It (i) selects the most discriminative CpGs per disorder class vs controls (|Δβ| ≥ `--min-delta`
and M-value t-test p < `--alpha`, top `--features-per-class`), (ii) builds class centroids,
and (iii) runs **nested leave-one-out cross-validation** (features re-selected each fold — no
leakage) reporting per-class sensitivity/specificity and overall accuracy. Outputs a bundle:
`reference_matrix.tsv` + `labels.tsv` (drop-in for `--episign-ref/--episign-labels`),
`centroids.tsv`, `signatures.json`, and `manifest.json` (params + CV metrics). Use `--no-validate`
to skip CV.

### `methylforge array` — Illumina EPIC/450K ingestion
Maps a methylation-array beta matrix into the same coordinate-keyed schema as sequencing data,
so array-based episignatures/imprinting use the identical downstream. The IDAT→beta step runs
in R (`bin/epic_process.R`, `modules/stage0_array.nf`, env `array.yml`, minfi funnorm + detP QC),
emitting a probe×sample beta matrix and a probe→coordinate manifest.
```bash
# probe x sample beta (+ manifest from minfi's EPIC annotation) -> coord matrix + per-sample methyl
methylforge array --beta epic_beta.tsv --manifest epic_manifest.tsv \
  --out-matrix epic_coord_matrix.tsv --out-methyl-dir epic_methyl/
# then, e.g.:  methylforge run epic_methyl/SAMPLE.methyl.tsv -s SAMPLE ...
```
Arrays carry no read depth, so a nominal coverage (`--cov`, default 30) is assigned to each beta
(flagged in QC). If the beta matrix is already coordinate-keyed, omit `--manifest`.

### `methylforge testdata` — synthetic fixture
```bash
methylforge testdata -o mf_demo
```

---

## 5. Nextflow usage

```bash
export NXF_SYNTAX_PARSER=v1     # the classic-DSL2 modules use Nextflow's v1 parser
nextflow run main.nf -profile <exec>,<pkg> --mode <mode> [options]
```
Profiles: exec ∈ `local|hpc_slurm`; pkg ∈ `conda|mamba|docker|singularity|apptainer`.

**Entry modes (`--mode`):**
- `sequencing` (default) — short/long/precalled methylation → per-sample analysis + reports, then an automatic batch QC dashboard. Needs `--samplesheet` (+ `--reference` for FASTQ/modBAM input).
- `array` — Illumina EPIC/450K. Either from raw IDATs (`--idat_dir` + `--array_sheet`, or `--array_sheet auto` to use a sample sheet inside the IDAT dir) → beta matrix via minfi (funnorm; noob for a single sample) with hg38 coordinates via sesame (`--array_genome hg38|hg19`); **or** from a pre-computed matrix (`--array_beta beta.tsv --array_manifest manifest.tsv`, skipping the R step). Then per-sample methylation → the same analysis + reports. The IDAT→beta step needs the R env (`env/array.yml`); validate it on the cluster with `bash test/test_array_hpc.sh` (uses minfi's bundled example IDATs — no patient data).
- `cohort_qc` — batch QC dashboard only, from `--samplesheet` (`sample,methylation`).

Key params: `--targets` (targeted BED), `--library bisulfite|emseq`, `--aligner bismark|bwameth`, `--episign_ref/--episign_labels`, `--run_dmr --dmr_matrix --dmr_groups`, `--genes_bed`, `--report_modes`, `--cohort_qc true|false` (batch dashboard), `--max_cpus/--max_memory`.

Every `sequencing`/`array` run also writes `${outdir}/cohort_qc.html` aggregating all samples. The `test` profile runs the analysis+report engine on the synthetic fixture (no aligners): `nextflow run main.nf -profile test,conda`.

---

## 6. Pipeline stages

| Stage | Module | Tool(s) | Purpose |
|---|---|---|---|
| 1 QC/trim | `stage1_qc_trim.nf` | FastQC, fastp/Trim Galore | read QC, adapter/quality trimming (bisulfite-aware) |
| 2 align (short) | `stage2_align_short.nf` | Bismark or bwa-meth + dedup | methylation-aware alignment |
| 2 align (long) | `stage2_align_long.nf` | minimap2 `-y` + modkit | align native-5mC reads, CpG pileup |
| 3 extract | `stage3_extract.nf` | MethylDackel | per-CpG methylation (+CHH/CHG for conversion QC) |
| 4 QC | `qc.py` | — | conversion, coverage, global methylation, sex, gates |
| 5 imprinting | `imprinting.py` | — | imprinted-DMR loss/gain → disorder, MLID |
| 6 DMR | `dmr.py` | — | DMP (Mann-Whitney+BH) → merged DMRs |
| 7 episignature | `episignature.py` | (scikit-learn opt.) | reference-based classification |
| 8 annotate | `annotate.py` | (annotatr in prod) | nearest-gene / context for DMRs |
| 9 report | `report.py` | Jinja2, matplotlib | clinical + research HTML, results.json |

---

## 7. Report engine

Two styles from one result set:
- **clinical** → `clinical_gbcf.html.j2`: GBCF/Noguchi-branded, A4 print CSS, logos in the header, a **result summary** table, imprinted-locus plot+table, episignature, QC, an interpretation/limitations note, and Analysed-by/Authorised-by signature lines. Findings are coloured by the **ACMG/AMP severity convention**: Pathogenic (red), Likely pathogenic (orange→red), Uncertain/VUS (amber), Benign/normal (green). Severity mapping: an imprinted-locus loss/gain with a mapped disorder → *Pathogenic*; aberrant without a mapped disorder → *Likely pathogenic*; a positive episignature → *Pathogenic*; inconclusive episignature → *VUS*; normal/control-like → *Benign*.
- **research** → `report.html.j2`: modern theme-aware layout with the full DMR table + volcano.

Branding assets are bundled in `methylforge/resources/assets/` (`gbcf_lab_logo.png`, `noguchi_logo.png`) and embedded as base64, so reports are self-contained on any machine.

---

## 8. FamilyForge interoperability

Every `run` emits `SAMPLE.familyforge.json` — a `candidates.json`-compatible object whose `rows` use the same base keys as FamilyForge candidates (`gene, model, tier, tier_reason, confidence, segregation, zygosity, flags`) plus epigenetic extras (`finding_type, region, call, disorder, metric`). Tiering mirrors FamilyForge A/B/C: imprinting defect or positive episignature → **A**; DMR over a disease gene → **B**; research DMRs → **C**. A companion tolerant `SAMPLE.epi_findings.tsv` (gene-keyed) can be ingested like CallForge's `cnv_calls.tsv`. This lets epigenetic findings flow into FamilyForge's interpretation and clinical-report generation alongside sequence variants.

---

## 9. QC gates & interpretation

Defaults (override via the API / config): conversion efficiency ≥ 0.98, mean CpG coverage ≥ 10×, ≥ 1000 callable CpGs (≥ 5× each), pUC19 methylated-control retention ≥ 0.95. A failing sample is flagged **REVIEW**, never silently dropped.

**Bisulfite/EM conversion controls (spike-ins).** Conversion is estimated, in order of preference, from (1) **spike-in controls auto-detected in the aligned data** — *lambda phage* (unmethylated → the conversion rate) and *CpG-methylated pUC19* (→ the over-conversion / false-negative rate, gated separately); (2) a provided value; or (3) non-CpG (CHH) methylation. The spike-in FASTAs are downloaded to `_local/refs/spikein_controls.fa`; append them to your genome once so spike-in reads align and get methylation-called:
```bash
bin/prepare_reference.sh GRCh38.fa            # -> GRCh38.withspikein.fa, then index it
```
MethylForge auto-detects the lambda (`NC_001416`) and pUC19 (`L09137`) contigs in the calls — no extra flags. Contig names are configurable via the `spikein=` argument to `methylation_qc`. Spike-in contigs are excluded from the biological CpG metrics.

---

## 10. Best practices & standards

- **Conversion controls:** include an unmethylated spike-in (lambda) and, where possible, a fully-methylated control; report conversion efficiency on every run.
- **Coverage:** imprinting/targeted calls need adequate per-locus depth; genome-wide WGBS is data-heavy — plan storage/compute accordingly.
- **Confirmation:** confirm clinically significant imprinting defects with an orthogonal assay (MS-MLPA is the clinical workhorse) before diagnostic reporting.
- **Ancestry-matched references:** episignature classifiers and methylation reference ranges are overwhelmingly European-ancestry; build/validate local controls — this is both a rigour requirement and a research opportunity.
- **Provenance & reproducibility:** pinned per-stage conda envs / containers; every run records provenance in `results.json`.

---

## 11. Limitations & validation

- The bundled imprinted-DMR panel coordinates/thresholds and the episignature reference are **starting references** and must be locally validated.
- MethylForge reports are **research-grade**; they are not stand-alone diagnoses.
- Sequence-based methylation does not, by itself, resolve parent-of-origin/UPD without additional data; long-read phasing helps.

---

## 12. Extending

- **Add an imprinting locus:** append a row to `resources/imprinted_dmrs.tsv`.
- **Add a disorder episignature:** add reference samples + a label to the episignature reference/labels; the classifier picks it up automatically.
- **Swap aligners:** `--aligner bismark|bwameth`; add EM-seq via `--library emseq`.

---

## 13. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| conversion `n/a` | no CHH sites / spike-in | run MethylDackel with `--CHH`, or pass a spike-in value |
| everything "no_call" at imprinted loci | low coverage in those regions | check targeted BED / depth; imprinting needs ≥3 covered CpGs/locus |
| episignature always "inconclusive" | reference too small / features absent in sample | enlarge reference; ensure feature CpGs are covered |
| Apple-Silicon env solve fails | tool lacks arm64 build | `CONDA_SUBDIR=osx-64 mamba env create -f env/<x>.yml` |
