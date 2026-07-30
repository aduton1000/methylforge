"""MethylForge orchestration — assemble the analysis stages for one sample (+ optional
cohort inputs) and emit clinical + research reports and a results.json."""
from __future__ import annotations
import os
import json
import datetime
import numpy as np
import pandas as pd
from . import io as eio, qc as eqc, imprinting as eimp, dmr as edmr, \
    episignature as eepi, annotate as eann, report as erep


def _query_from_meth(meth: pd.DataFrame, features) -> pd.Series:
    lut = {f"{c}:{s}": b for c, s, b in zip(meth.chrom, meth.start, meth.beta)}
    return pd.Series({f: lut.get(f, np.nan) for f in features}, dtype=float)


def run(meth_path: str, sample: str, outdir: str, *,
        family: str | None = None,
        assay: str = "methylation", genome: str = "GRCh38",
        modes=("clinical", "research"),
        panel_path: str | None = None,
        episign_ref: str | None = None, episign_labels: str | None = None,
        dmr_matrix: str | None = None, dmr_groups: str | None = None,
        genes_bed: str | None = None, gates: dict | None = None,
        lambda_meth: float | None = None, excel: bool = True) -> dict:
    os.makedirs(outdir, exist_ok=True)
    meth = eio.read_methylation(meth_path)

    # Stage 4 — QC
    qc = eqc.methylation_qc(meth, gates=gates, lambda_meth=lambda_meth)

    # Stage 5 — imprinting
    panel = eio.read_dmr_panel(panel_path)
    imprinting = eimp.analyse_imprinting(meth, panel)
    imprint_summary = eimp.summarise_imprinting(imprinting)

    results = dict(sample=sample, family=family or sample, assay=assay, genome=genome, meth=meth, qc=qc,
                   imprinting=imprinting, imprint_summary=imprint_summary,
                   provenance=f"methylforge run · {datetime.datetime.now().isoformat(timespec='seconds')} · input={os.path.basename(meth_path)}")

    # Stage 7 — episignature (optional)
    if episign_ref and episign_labels:
        ref = pd.read_csv(episign_ref, sep="\t", index_col=0)
        lab = pd.read_csv(episign_labels, sep="\t")
        labels = dict(zip(lab.iloc[:, 0], lab.iloc[:, 1]))
        centroids = eepi.build_centroids(ref, labels)
        query = _query_from_meth(meth, ref.index)
        results["episign"] = eepi.classify(query, centroids, ref, labels)

    # Stage 6 — DMR (optional, cohort)
    if dmr_matrix and dmr_groups:
        mat = pd.read_csv(dmr_matrix, sep="\t", index_col=0)
        gdf = pd.read_csv(dmr_groups, sep="\t")
        groups = dict(zip(gdf.iloc[:, 0], gdf.iloc[:, 1]))
        dmps = edmr.call_dmps(mat, groups)
        dmrs = edmr.merge_dmrs(dmps)
        genes = eann.load_genes_bed(genes_bed) if genes_bed else None
        dmrs = eann.annotate_regions(dmrs, genes)
        results.update(dmps=dmps, dmrs=dmrs, dmr_summary=edmr.summarise_dmr(dmps, dmrs))

    # Stage 9 — reports
    outputs = {}
    for mode in modes:
        out_html = os.path.join(outdir, f"{sample}.{mode}.report.html")
        erep.render_report(results, out_html, mode=mode)
        outputs[mode] = out_html
    erep.write_results_json(results, os.path.join(outdir, f"{sample}.results.json"))

    # Stage 9b — Excel findings workbook (Findings · Imprinting loci · QC)
    if excel:
        from . import export_excel as exxl
        xlsx = os.path.join(outdir, f"{sample}.findings.xlsx")
        exxl.write_sample_findings_xlsx(results, xlsx)
        outputs["findings_xlsx"] = xlsx

    results["outputs"] = outputs
    return results
