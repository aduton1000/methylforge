#!/usr/bin/env bash
# Append bisulfite-conversion spike-in controls to the reference genome so that
# spike-in reads align and MethylForge's QC can measure conversion efficiency
# (lambda) and over-conversion (CpG-methylated pUC19) from the aligned data.
# Usage: bin/prepare_reference.sh <genome.fa> [<spikein.fa>] [<out.fa>]
set -euo pipefail
GENOME=$1
SPIKEIN=${2:-"$(dirname "$0")/../_local/refs/spikein_controls.fa"}
OUT=${3:-"${GENOME%.fa*}.withspikein.fa"}
[ -f "$GENOME" ] || { echo "genome not found: $GENOME" >&2; exit 1; }
[ -f "$SPIKEIN" ] || { echo "spike-in FASTA not found: $SPIKEIN (run the DB download)" >&2; exit 1; }
cat "$GENOME" "$SPIKEIN" > "$OUT"
echo "wrote $OUT"
echo "Next: index it (samtools faidx; bismark_genome_preparation <dir>; or bwameth.py index),"
echo "then MethylForge QC auto-detects the lambda (NC_001416) / pUC19 (L09137) contigs."
