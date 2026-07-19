#!/usr/bin/env bash
# End-to-end demo of the MethylForge analysis + report engine on synthetic data.
set -euo pipefail
OUT=${1:-mf_demo}
echo "[1/5] generating synthetic test data ..."
methylforge testdata -o "$OUT"
D="$OUT/data"
echo "[2/5] PWS case (expect Prader-Willi imprinting call) ..."
methylforge run "$D/pws_case.methyl.tsv" -s pws_case --family FAM01 \
  --episign-ref "$D/episig_reference.tsv" --episign-labels "$D/episig_labels.tsv" -o "$OUT/pws"
echo "[3/5] Sotos case (expect positive episignature + DMR) ..."
methylforge run "$D/sotos_case.methyl.tsv" -s sotos_case --family FAM02 \
  --episign-ref "$D/episig_reference.tsv" --episign-labels "$D/episig_labels.tsv" \
  --dmr-matrix "$D/dmr_matrix.tsv" --dmr-groups "$D/dmr_groups.tsv" --genes-bed "$D/genes.bed" -o "$OUT/sotos"
echo "[4/5] normal control (expect all normal) ..."
methylforge run "$D/normal_control.methyl.tsv" -s normal_control \
  --episign-ref "$D/episig_reference.tsv" --episign-labels "$D/episig_labels.tsv" -o "$OUT/normal"
echo "[5/5] done. Reports:"
find "$OUT" -name "*.report.html" | sort
