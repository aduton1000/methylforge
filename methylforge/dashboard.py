"""Cohort QC dashboard — a MultiQC-style batch overview.

Aggregates per-sample methylation QC across a whole run into one self-contained
HTML: a batch summary, overview plots (conversion, coverage, callable CpGs, global
methylation, spike-in controls, pass/fail), and a sortable per-sample table with
pass/fail badges and per-metric colouring. Also writes a summary TSV.
"""
from __future__ import annotations
import os
import json
import datetime
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
from . import io as eio, qc as eqc, plots, __version__

_TPL_DIR = os.path.join(os.path.dirname(__file__), "report_templates")


def _row_from_qc(name, m):
    spk = m.get("spikein", {}) or {}
    return dict(
        sample=name, qc_pass=bool(m["qc_pass"]),
        conversion=m["conversion_efficiency"], conversion_source=m["conversion_source"],
        mean_cov=m["mean_cpg_coverage"], median_cov=m["median_cpg_coverage"],
        callable=m["n_cpg_callable"], n_cpg=m["n_cpg"],
        global_cpg_meth=m["global_cpg_methylation"], global_chh=m.get("global_chh_methylation"),
        lambda_conversion=spk.get("lambda_conversion"),
        puc19_methylation=spk.get("puc19_methylation"),
        puc19_overconversion=spk.get("puc19_overconversion"),
        inferred_sex=(m["inferred_sex"].split()[0] if m.get("inferred_sex") else "unknown"),
        checks=m.get("checks", {}))


def collect_qc(sample_paths: dict, gates: dict | None = None) -> pd.DataFrame:
    rows = []
    for name, path in sample_paths.items():
        m = eqc.methylation_qc(eio.read_methylation(path), gates=gates)
        rows.append(_row_from_qc(name, m))
    return pd.DataFrame(rows)


def collect_from_results(result_paths) -> pd.DataFrame:
    """Aggregate QC already computed in per-sample results.json files."""
    rows = []
    for p in result_paths:
        d = json.load(open(p))
        rows.append(_row_from_qc(d.get("sample", os.path.basename(p)), d["qc"]))
    return pd.DataFrame(rows)


def render_dashboard(df: pd.DataFrame, out_html: str, title="MethylForge cohort QC",
                     gates: dict | None = None, date: str | None = None) -> str:
    g = {**eqc.DEFAULT_GATES, **(gates or {})}
    env = Environment(loader=FileSystemLoader(_TPL_DIR),
                      autoescape=select_autoescape(["html", "xml"]))
    tpl = env.get_template("cohort_qc.html.j2")
    df = df.sort_values("sample").reset_index(drop=True)

    def fail_cells(r):
        c = r["checks"] or {}
        return dict(conversion=not c.get("conversion", True),
                    coverage=not c.get("coverage", True),
                    callable=not c.get("callable", True),
                    overconversion=not c.get("overconversion", True))

    rows = []
    for _, r in df.iterrows():
        rows.append({**r.to_dict(), "fail": fail_cells(r)})

    html = tpl.render(
        title=title, version=__version__,
        date=date or datetime.date.today().isoformat(),
        n=len(df), n_pass=int(df["qc_pass"].sum()), n_review=int((~df["qc_pass"]).sum()),
        gates=g, rows=rows, figure=plots.cohort_qc_figure(df, g),
    )
    with open(out_html, "w") as fh:
        fh.write(html)
    return out_html


def write_summary(df: pd.DataFrame, path: str) -> str:
    cols = ["sample", "qc_pass", "conversion", "conversion_source", "mean_cov", "callable",
            "global_cpg_meth", "lambda_conversion", "puc19_overconversion", "inferred_sex"]
    df[[c for c in cols if c in df.columns]].to_csv(path, sep="\t", index=False)
    return path
