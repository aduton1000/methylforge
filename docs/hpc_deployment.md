# MethylForge on a shared Ubuntu HPC

This guide installs MethylForge once on the lab's Ubuntu Slurm cluster so any member can run it with a single command.

## 1. Prerequisites (per cluster)
- Nextflow ≥ 23.10 (Java 17+): `curl -s https://get.nextflow.io | bash` (or an `lmod` module).
- Micromamba/Mamba, **or** SingularityCE/Apptainer (recommended for multi-user reproducibility).
- Shared, backed-up storage for reference genomes, episignature references, and results.

## 2. Clone to a shared location
```bash
cd /opt/lab/pipelines            # a path all users can read
git clone https://github.com/aduton1000/methylforge.git
cd methylforge
```

## 3. Environments — choose ONE strategy

### 3a. Conda/mamba (native on linux-64)
```bash
mamba env create -f env/methylforge.yml   # core engine
conda activate methylforge && pip install -e .
# heavy-aligner envs build on demand during a run (or pre-create):
mamba env create -f env/align_short.yml
mamba env create -f env/align_long.yml
mamba env create -f env/extract.yml
```
Run with `-profile hpc_slurm,mamba`.

### 3b. Singularity (recommended for reproducibility)
Build/pull one image per env and run with `-profile hpc_slurm,singularity`. Nextflow pulls the per-process `conda` directives into Singularity automatically when the `singularity` profile is active and images are configured; alternatively bake images with `apptainer build`. Point Nextflow at a shared cache:
```bash
export NXF_SINGULARITY_CACHEDIR=/opt/lab/singularity_cache
nextflow run main.nf -profile hpc_slurm,singularity --samplesheet samples.csv --reference $REF
```

## 4. Slurm profile
`nextflow.config` already defines `hpc_slurm` (executor=slurm, queue chosen by task time, queueSize=100). Adjust `process.queue`, `--max_cpus`, `--max_memory`, and per-label resources in `conf/base.config` to your partitions.

## 5. One-command wrapper for the lab
Create `/opt/lab/pipelines/methylforge/bin/methylforge-run`:
```bash
#!/usr/bin/env bash
set -euo pipefail
module load nextflow singularity 2>/dev/null || true
export NXF_SINGULARITY_CACHEDIR=/opt/lab/singularity_cache
REPO=/opt/lab/pipelines/methylforge
REF=/opt/lab/refs/GRCh38.full.fa
nextflow run "$REPO/main.nf" -profile hpc_slurm,singularity \
  --reference "$REF" \
  --episign_ref /opt/lab/refs/episig_reference.tsv \
  --episign_labels /opt/lab/refs/episig_labels.tsv \
  "$@"
```
`chmod +x` it and add to PATH (or ship an `lmod` module). A user then runs:
```bash
methylforge-run --samplesheet my_samples.csv --outdir /scratch/$USER/methyl_run
```

## 6. Reference data layout (shared)
```
/opt/lab/refs/
  GRCh38.full.fa (+ .fai; Bismark genome prep dir; bwa-meth index)
  episig_reference.tsv   episig_labels.tsv     # locally-validated, ancestry-matched
  imprinted_dmrs.tsv                            # curated + locally validated
  genes.bed
```
Build the Bismark index once (`bismark_genome_preparation <dir>`); MethylForge reuses it.

## 7. Storage & compute notes
- WGBS is data-heavy; give runs generous `/scratch` and clean Nextflow `work/` between large runs (`nextflow clean -f`).
- Long-read basecalling with modified-base models (ONT `dorado`) is **GPU** and done upstream; MethylForge consumes the resulting modBAM.
- Keep `NXF_SINGULARITY_CACHEDIR` and the conda pkgs cache on shared, fast storage to avoid per-user rebuilds.

## 8. Validation on the cluster
```bash
nextflow run main.nf -profile test,singularity     # synthetic fixture, analysis+report only
python tests/selftest.py                            # asserts expected calls
```
Both should pass before running patient data.
