"""Stage 4 — methylation quality control.

Computes the metrics that gate a methylation experiment, per international
best-practice (bisulfite conversion efficiency, coverage, global methylation
levels by context, and an inferred sex), and applies pass/fail thresholds.
Conversion efficiency is estimated from non-CpG (CHH) methylation when present
(1 - mean CHH beta) or from a spike-in control if supplied.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

DEFAULT_GATES = dict(
    min_conversion=0.98,      # bisulfite/EM conversion efficiency
    min_mean_cov=10.0,        # mean CpG coverage
    min_cpg_callable=1000,    # CpGs at >= min_site_cov
    min_site_cov=5,
)


def methylation_qc(meth: pd.DataFrame, gates: dict | None = None,
                   lambda_meth: float | None = None) -> dict:
    g = {**DEFAULT_GATES, **(gates or {})}
    cpg = meth[meth.context == "CpG"]
    chh = meth[meth.context == "CHH"]
    covered = cpg[cpg['cov'] > 0]
    callable_sites = int((cpg['cov'] >= g["min_site_cov"]).sum())

    if lambda_meth is not None:
        conversion = 1.0 - float(lambda_meth)
        conv_source = "spike-in control"
    elif len(chh) > 0:
        conversion = 1.0 - float(np.nansum(chh.meth) / max(np.nansum(chh['cov']), 1))
        conv_source = "non-CpG (CHH) methylation"
    else:
        conversion = np.nan
        conv_source = "not estimable (no CHH / spike-in)"

    metrics = {
        "n_cytosines": int(len(meth)),
        "n_cpg": int(len(cpg)),
        "n_cpg_covered": int(len(covered)),
        "n_cpg_callable": callable_sites,
        "mean_cpg_coverage": round(float(covered['cov'].mean()) if len(covered) else 0.0, 2),
        "median_cpg_coverage": float(covered['cov'].median()) if len(covered) else 0.0,
        "global_cpg_methylation": round(float(np.nansum(cpg.meth) / max(np.nansum(cpg['cov']), 1)), 4),
        "global_chg_methylation": round(float(np.nansum(meth[meth.context == 'CHG'].meth) /
                                              max(np.nansum(meth[meth.context == 'CHG']['cov']), 1)), 4)
        if (meth.context == "CHG").any() else None,
        "global_chh_methylation": round(float(np.nansum(chh.meth) / max(np.nansum(chh['cov']), 1)), 4)
        if len(chh) else None,
        "conversion_efficiency": round(conversion, 4) if conversion == conversion else None,
        "conversion_source": conv_source,
    }
    metrics["inferred_sex"] = _infer_sex(cpg)

    checks = {
        "conversion": (metrics["conversion_efficiency"] is None) or
                      (metrics["conversion_efficiency"] >= g["min_conversion"]),
        "coverage": metrics["mean_cpg_coverage"] >= g["min_mean_cov"],
        "callable": metrics["n_cpg_callable"] >= g["min_cpg_callable"],
    }
    metrics["checks"] = checks
    metrics["gates"] = g
    metrics["qc_pass"] = all(checks.values())
    return metrics


def _infer_sex(cpg: pd.DataFrame) -> str:
    """Rough sex inference from chrY coverage relative to autosomes."""
    chrom = cpg.chrom.str.replace("chr", "", regex=False)
    auto = cpg[chrom.isin([str(i) for i in range(1, 23)])]
    y = cpg[chrom == "Y"]
    if len(auto) == 0:
        return "unknown"
    auto_cov = auto['cov'].mean()
    y_cov = y['cov'].mean() if len(y) else 0.0
    ratio = (y_cov / auto_cov) if auto_cov else 0.0
    if ratio >= 0.15:
        return f"male (chrY/auto coverage ratio {ratio:.2f})"
    if ratio < 0.05:
        return f"female (chrY/auto coverage ratio {ratio:.2f})"
    return f"indeterminate (chrY/auto coverage ratio {ratio:.2f})"
