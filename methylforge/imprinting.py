"""Stage 5 — imprinted-DMR analysis.

Interrogates methylation at a curated panel of clinically important imprinting
control regions and calls loss/gain of methylation against per-locus normal
ranges, mapping the direction of the defect to the associated disorder
(e.g., SNRPN hypermethylation -> Prader-Willi; H19/IGF2 hypomethylation ->
Silver-Russell). Also flags multi-locus imprinting disturbance (MLID).

NOTE: panel coordinates and thresholds are a curated starting reference and MUST
be locally validated (ideally against ancestry-matched controls) before any
diagnostic use.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from . import io as eio

# call thresholds around the expected ~0.5 (one methylated allele)
HYPER = 0.65
HYPO = 0.35
MIN_SITES = 3     # minimum covered CpGs in the region to make a call


def _region_beta(meth: pd.DataFrame, chrom, start, end) -> tuple[float, int]:
    sub = meth[(meth.chrom == chrom) & (meth.start >= start) & (meth.end <= end) &
               (meth.context == "CpG") & (meth['cov'] > 0)]
    if len(sub) == 0:
        return (np.nan, 0)
    # coverage-weighted mean methylation
    beta = float(np.nansum(sub.meth) / max(np.nansum(sub['cov']), 1))
    return (beta, int(len(sub)))


def analyse_imprinting(meth: pd.DataFrame, panel: pd.DataFrame | None = None,
                       hyper=HYPER, hypo=HYPO) -> pd.DataFrame:
    if panel is None:
        panel = eio.read_dmr_panel()
    rows = []
    for _, r in panel.iterrows():
        beta, n = _region_beta(meth, r.chrom, int(r.start), int(r.end))
        call, disorder = "no_call", ""
        if n >= MIN_SITES and beta == beta:
            delta = beta - float(r.normal_meth)
            if beta >= hyper:
                call = "hypermethylation"
                disorder = str(r.hyper_disorder) if pd.notna(r.hyper_disorder) else ""
            elif beta <= hypo:
                call = "hypomethylation"
                disorder = str(r.hypo_disorder) if pd.notna(r.hypo_disorder) else ""
            else:
                call = "normal"
        else:
            delta = np.nan
        rows.append(dict(locus=r.locus, gene=r.gene, chrom=r.chrom,
                         start=int(r.start), end=int(r.end),
                         normal_meth=float(r.normal_meth),
                         observed_beta=round(beta, 4) if beta == beta else None,
                         n_cpg=n, delta=round(delta, 4) if delta == delta else None,
                         call=call, disorder=disorder, notes=r.notes))
    out = pd.DataFrame(rows)
    return out


def summarise_imprinting(res: pd.DataFrame) -> dict:
    abnormal = res[res.call.isin(["hypermethylation", "hypomethylation"])]
    disorders = sorted({d for d in abnormal.disorder if d})
    mlid = len(abnormal) >= 2
    if len(abnormal) == 0:
        headline = "No abnormal methylation detected at the imprinted loci examined."
        severity = "normal"
    elif mlid:
        headline = (f"Multi-locus imprinting disturbance (MLID): {len(abnormal)} loci "
                    f"abnormal — {', '.join(abnormal.locus)}.")
        severity = "abnormal"
    else:
        a = abnormal.iloc[0]
        headline = (f"{a.call} at {a.locus} ({a.gene}) — consistent with "
                    f"{a.disorder or 'an imprinting disorder'}.")
        severity = "abnormal"
    return dict(headline=headline, severity=severity, n_abnormal=int(len(abnormal)),
                disorders=disorders, mlid=bool(mlid))
