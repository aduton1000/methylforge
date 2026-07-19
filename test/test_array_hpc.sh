#!/usr/bin/env bash
# HPC self-test of the array IDAT -> beta -> per-sample methylation chain, using
# minfi's bundled example 450K IDATs (no patient data needed).
#   conda activate methylforge-array   # env/array.yml  (must also have `methylforge` on PATH)
#   bash test/test_array_hpc.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-mf_array_test}; mkdir -p "$OUT"; cd "$OUT"
echo "[1/4] locating minfiData example IDATs ..."
IDAT=$(Rscript -e 'cat(system.file("extdata", package="minfiData"))')
[ -n "$IDAT" ] && [ -d "$IDAT" ] || { echo "minfiData not installed — create env/array.yml"; exit 1; }
echo "      $IDAT"
echo "[2/4] IDAT -> beta + manifest (bin/epic_process.R, auto sample sheet, hg38) ..."
CHAIN="$HERE/_local/refs/hg19ToHg38.over.chain"   # deterministic hg19->hg38 fallback (optional)
Rscript "$HERE/bin/epic_process.R" "$IDAT" auto epic_beta.tsv epic_qc.txt hg38 "$([ -f "$CHAIN" ] && echo "$CHAIN")"
echo "----- epic_qc.txt -----"; cat epic_qc.txt
echo "[3/4] beta + manifest -> per-sample methylation (methylforge array) ..."
methylforge array --beta epic_beta.tsv --manifest epic_manifest.tsv \
  --out-matrix coord_matrix.tsv --out-methyl-dir methyl
echo "[4/4] QC one converted sample ..."
S=$(ls methyl/*.methyl.tsv | head -1)
methylforge qc "$S" | head -18
echo "OK — array IDAT -> beta -> methylation chain works on this machine."
