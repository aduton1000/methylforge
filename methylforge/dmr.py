"""Stage 6 — differential methylation (DMP / DMR) calling.

Case-vs-control (or group-vs-group) differential methylation on a per-CpG
methylation matrix. Per-CpG differentially methylated positions (DMPs) are
detected with a Mann-Whitney U test on beta values plus a mean methylation
difference threshold; adjacent DMPs are merged into differentially methylated
regions (DMRs). Multiple-testing correction is Benjamini-Hochberg.

Input matrix: rows = CpG sites (index like 'chrom:pos'), columns = samples,
values = beta (0..1). A grouping maps samples -> 'case' / 'control'.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def call_dmps(matrix: pd.DataFrame, groups: dict, min_delta=0.10,
              fdr=0.05) -> pd.DataFrame:
    cases = [s for s, g in groups.items() if g == "case" and s in matrix.columns]
    ctrls = [s for s, g in groups.items() if g == "control" and s in matrix.columns]
    if not cases or not ctrls:
        raise ValueError("need >=1 case and >=1 control sample present in the matrix")
    C = matrix[cases].to_numpy(float)
    K = matrix[ctrls].to_numpy(float)
    delta = np.nanmean(C, axis=1) - np.nanmean(K, axis=1)
    pvals = np.ones(len(matrix))
    for i in range(len(matrix)):
        a, b = C[i][~np.isnan(C[i])], K[i][~np.isnan(K[i])]
        if len(a) >= 1 and len(b) >= 1 and (len(a) + len(b)) >= 3:
            try:
                if len(a) >= 2 and len(b) >= 2:
                    pvals[i] = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
                else:
                    pvals[i] = stats.ttest_ind(a, b, equal_var=False).pvalue
            except Exception:
                pvals[i] = 1.0
    q = _bh_fdr(pvals)
    out = pd.DataFrame({
        "site": matrix.index, "delta": np.round(delta, 4),
        "case_mean": np.round(np.nanmean(C, axis=1), 4),
        "control_mean": np.round(np.nanmean(K, axis=1), 4),
        "pvalue": pvals, "qvalue": np.round(q, 5)})
    out["chrom"] = out.site.str.split(":").str[0]
    out["pos"] = out.site.str.split(":").str[1].astype(int)
    out["dmp"] = (np.abs(out.delta) >= min_delta) & (out.qvalue <= fdr)
    out["direction"] = np.where(out.delta > 0, "hyper", "hypo")
    return out.sort_values(["chrom", "pos"]).reset_index(drop=True)


def merge_dmrs(dmps: pd.DataFrame, max_gap=1000, min_cpg=3) -> pd.DataFrame:
    """Merge significant, same-direction DMPs within max_gap bp into DMRs."""
    sig = dmps[dmps.dmp].copy()
    regions = []
    for (chrom, direction), grp in sig.groupby(["chrom", "direction"]):
        grp = grp.sort_values("pos")
        start = prev = None
        buf = []
        for _, row in grp.iterrows():
            if start is None:
                start = row.pos; buf = [row]
            elif row.pos - prev <= max_gap:
                buf.append(row)
            else:
                regions.append(_finalise_dmr(chrom, direction, buf))
                start = row.pos; buf = [row]
            prev = row.pos
        if buf:
            regions.append(_finalise_dmr(chrom, direction, buf))
    dmr = pd.DataFrame(regions)
    if len(dmr):
        dmr = dmr[dmr.n_cpg >= min_cpg].sort_values("mean_abs_delta", ascending=False)
        dmr = dmr.reset_index(drop=True)
    return dmr


def _finalise_dmr(chrom, direction, buf):
    rows = pd.DataFrame(buf)
    return dict(chrom=chrom, start=int(rows.pos.min()), end=int(rows.pos.max()),
                width=int(rows.pos.max() - rows.pos.min() + 1),
                n_cpg=int(len(rows)), direction=direction,
                mean_delta=round(float(rows.delta.mean()), 4),
                mean_abs_delta=round(float(rows.delta.abs().mean()), 4),
                min_qvalue=float(rows.qvalue.min()))


def summarise_dmr(dmps: pd.DataFrame, dmrs: pd.DataFrame) -> dict:
    n_dmp = int(dmps.dmp.sum())
    n_dmr = int(len(dmrs))
    hyper = int((dmrs.direction == "hyper").sum()) if n_dmr else 0
    hypo = int((dmrs.direction == "hypo").sum()) if n_dmr else 0
    return dict(n_tested=int(len(dmps)), n_dmp=n_dmp, n_dmr=n_dmr,
                n_hyper=hyper, n_hypo=hypo,
                headline=f"{n_dmr} DMRs ({hyper} hyper-, {hypo} hypomethylated) "
                         f"from {n_dmp} significant CpGs.")
