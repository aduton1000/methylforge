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


def _full_query(meth: pd.DataFrame) -> pd.Series:
    q = pd.Series({f"{c}:{s}": b for c, s, b in zip(meth.chrom, meth.start, meth.beta)}, dtype=float)
    return q[~q.index.duplicated(keep="first")]


def run(meth_path: str, sample: str, outdir: str, *,
        family: str | None = None,
        assay: str = "methylation", genome: str = "GRCh38",
        modes=("clinical", "research"),
        panel_path: str | None = None,
        episign_ref: str | None = None, episign_labels: str | None = None,
        dmr_matrix: str | None = None, dmr_groups: str | None = None,
        genes_bed: str | None = None, gates: dict | None = None,
        lambda_meth: float | None = None, excel: bool = True,
        gates_profile: str = "standard",
        cell_reference: str | None = None, decon_method: str = "CP",
        clocks: dict | None = None, chrono_age: float | None = None,
        dmr_covariates: str | None = None, dmr_pairs: str | None = None,
        dmr_autosomes: bool = False,
        asm_hp1: str | None = None, asm_hp2: str | None = None) -> dict:
    os.makedirs(outdir, exist_ok=True)
    meth = eio.read_methylation(meth_path)

    # Stage 4 — QC (gates_profile='low-pass' relaxes coverage gates for shallow data)
    if gates is None and gates_profile != "standard":
        gates = eqc.gates_for(gates_profile)
    qc = eqc.methylation_qc(meth, gates=gates, lambda_meth=lambda_meth)

    # Stage 5 — imprinting
    panel = eio.read_dmr_panel(panel_path)
    imprinting = eimp.analyse_imprinting(meth, panel)
    imprint_summary = eimp.summarise_imprinting(imprinting)

    results = dict(sample=sample, family=family or sample, assay=assay, genome=genome, meth=meth, qc=qc,
                   imprinting=imprinting, imprint_summary=imprint_summary,
                   provenance=f"methylforge run · {datetime.datetime.now().isoformat(timespec='seconds')} · input={os.path.basename(meth_path)}")

    # Stage 4b — cell-type deconvolution (optional; needs a validated reference)
    if cell_reference:
        from . import deconvolve as edc
        try:
            dref = edc.load_reference(cell_reference)
            dres = edc.deconvolve(_full_query(meth), dref, method=decon_method)
            results["deconvolution"] = {**dres, "summary": edc.summarise(dres)}
        except Exception as e:
            results["deconvolution"] = {"error": f"{type(e).__name__}: {e}"}

    # Stage 4c — epigenetic clocks (optional; clocks = {name: coefficient_path})
    if clocks:
        from . import clocks as eclk
        q = _full_query(meth)
        rows = []
        for name, path in clocks.items():
            ck = eclk.load_clock(path)
            ck["meta"]["name"] = name
            r = eclk.predict_age(q, ck)
            if chrono_age is not None:
                r["delta_age"] = eclk.delta_age(r["dnam_age"], chrono_age)
            rows.append(r)
        results["clocks"] = rows

    # Stage 7 — episignature (optional)
    if episign_ref and episign_labels:
        ref = pd.read_csv(episign_ref, sep="\t", index_col=0)
        lab = pd.read_csv(episign_labels, sep="\t")
        labels = dict(zip(lab.iloc[:, 0], lab.iloc[:, 1]))
        centroids = eepi.build_centroids(ref, labels)
        query = _query_from_meth(meth, ref.index)
        results["episign"] = eepi.classify(query, centroids, ref, labels)

    # Stage 6 — DMR (optional, cohort); covariate-adjusted / paired when supplied
    if dmr_matrix and dmr_groups:
        mat = pd.read_csv(dmr_matrix, sep="\t", index_col=0)
        gdf = pd.read_csv(dmr_groups, sep="\t")
        groups = dict(zip(gdf.iloc[:, 0], gdf.iloc[:, 1]))
        if dmr_autosomes:
            mat = edmr.restrict_autosomes(mat)
        covariates = pd.read_csv(dmr_covariates, sep=None, engine="python", index_col=0) \
            if dmr_covariates else None
        pair = None
        if dmr_pairs:
            pdf = pd.read_csv(dmr_pairs, sep=None, engine="python")
            pair = dict(zip(pdf.iloc[:, 0], pdf.iloc[:, 1]))
        if covariates is not None or pair is not None:
            dmps = edmr.call_dmps_glm(mat, groups, covariates=covariates, pair=pair)
        else:
            dmps = edmr.call_dmps(mat, groups)
        dmrs = edmr.merge_dmrs(dmps)
        genes = eann.load_genes_bed(genes_bed) if genes_bed else None
        dmrs = eann.annotate_regions(dmrs, genes)
        results.update(dmps=dmps, dmrs=dmrs, dmr_summary=edmr.summarise_dmr(dmps, dmrs))

    # Stage 6b — allele-specific methylation (optional; two haplotype methylation files)
    if asm_hp1 and asm_hp2:
        from . import asm as easm
        hp1, hp2 = easm.load_pair(asm_hp1, asm_hp2)
        asm_sites = easm.asm_from_haplotypes(hp1, hp2)
        asm_regions = easm.merge_asm_regions(asm_sites)
        results.update(asm=asm_regions, asm_summary=easm.summarise_asm(asm_sites, asm_regions))

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
