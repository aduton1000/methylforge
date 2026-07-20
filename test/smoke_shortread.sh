#!/usr/bin/env bash
# MethylForge short-read bisulfite SMOKE TEST — confirms the FASTQ -> methylation
# toolchain and a prepared bismark index work end-to-end on real tools:
#   simulate (Sherman) -> bismark align -> deduplicate -> MethylDackel -> methylforge run
#
# Prereqs (see docs/hpc_deployment.md / MANUAL install):
#   - conda envs: methylforge-align-short (bismark+bowtie2+samtools), methylforge-extract
#     (MethylDackel+samtools), and `methylforge` on PATH (from the methylforge env)
#   - a bismark-prepared genome directory (run `bismark_genome_preparation <dir>` once)
#
# Usage: bash test/smoke_shortread.sh <bismark_genome_dir> [n_reads] [outdir]
# NOTE: this is a plumbing check — with a small simulated read set the QC will flag
# low coverage (that's expected); the point is that every stage runs and reports.
set -eo pipefail
GDIR=${1:?usage: smoke_shortread.sh <bismark-prepared genome dir> [n_reads] [outdir]}
NREADS=${2:-120000}
W=${3:-mf_smoke}
GENOME=$(ls "$GDIR"/*.fa "$GDIR"/*.fasta 2>/dev/null | head -1)
[ -n "$GENOME" ] || { echo "no .fa/.fasta in $GDIR"; exit 1; }
SHORT="conda run --no-capture-output -n methylforge-align-short"
EXTRACT="conda run --no-capture-output -n methylforge-extract"

rm -rf "$W"; mkdir -p "$W/simsrc"; cd "$W"
echo "[0] fetch Sherman (bisulfite read simulator)"
[ -s Sherman ] || curl -sSL https://raw.githubusercontent.com/FelixKrueger/Sherman/master/Sherman -o Sherman

echo "[1] simulation source = first genome contig + lambda spike-in"
CONTIG=$($EXTRACT samtools faidx "$GENOME" 2>/dev/null; head -1 "$GENOME.fai" | cut -f1)
$EXTRACT samtools faidx "$GENOME" "$CONTIG" > simsrc/contig.fa
$EXTRACT samtools faidx "$GENOME" NC_001416.1 > simsrc/lambda.fa 2>/dev/null || true

echo "[2] Sherman: $NREADS paired-end 100 bp bisulfite reads"
$SHORT perl "$PWD/Sherman" --genome_folder "$PWD/simsrc" --length 100 \
  --number_of_seqs "$NREADS" --paired_end --error_rate 0.1 --output_dir "$PWD" >/dev/null 2>&1

echo "[3] bismark alignment"
$SHORT bismark --genome "$GDIR" -1 simulated_1.fastq -2 simulated_2.fastq -p 4 --basename smoke -o "$PWD" \
  2>&1 | grep -iE "mapping efficiency|unique best hit|completed"

echo "[4] deduplicate + [5] MethylDackel extract"
$SHORT deduplicate_bismark -p --bam --output_dir "$PWD" smoke_pe.bam 2>&1 | grep -iE "leftover" || true
$EXTRACT samtools sort -@4 -o smoke.sorted.bam smoke_pe.deduplicated.bam
$EXTRACT samtools index smoke.sorted.bam
$EXTRACT MethylDackel extract --CHH --CHG "$GENOME" smoke.sorted.bam -o smoke >/dev/null 2>&1
echo "  CpG sites called: $(($(wc -l < smoke_CpG.bedGraph) - 1))"

echo "[6] methylforge run"
methylforge run smoke_CpG.bedGraph -s SMOKE_TEST -o "$PWD/report" 2>&1 | tail -3
echo "OK — short-read bisulfite path (FASTQ -> bismark -> MethylDackel -> methylforge) works."
