"""Stage 4 — methylation quality control.

Computes the metrics that gate a methylation experiment, per international
best-practice (bisulfite conversion efficiency, coverage, global methylation
levels by context, and an inferred sex), and applies pass/fail thresholds.
Conversion efficiency is estimated, in order of preference, from (1) bisulfite/EM
spike-in controls auto-detected in the aligned data — lambda phage (unmethylated ->
conversion rate) and CpG-methylated pUC19 (-> over-conversion rate); (2) a provided
spike-in value; or (3) non-CpG (CHH) methylation (1 - mean CHH beta).
"""
from __future__ import annotations
import numpy as np
import pandas as pd

DEFAULT_GATES = dict(
    min_conversion=0.98,      # bisulfite/EM conversion efficiency
    min_mean_cov=10.0,        # mean CpG coverage
    min_cpg_callable=1000,    # CpGs at >= min_site_cov
    min_site_cov=5,
    min_puc19_meth=0.95,      # methylated-control retention (1 - over-conversion)
)

# spike-in contig names auto-detected in the methylation calls (case-insensitive,
# prefix match) — covers the downloaded accessions and common aliases.
DEFAULT_SPIKEIN = dict(
    lambda_contigs=["lambda", "chrlambda", "chrl", "nc_001416", "j02459",
                    "enterobacteria_phage_lambda"],
    puc19_contigs=["puc19", "chrpuc19", "l09137", "m77789"],
)


def spikein_conversion(meth: pd.DataFrame, spikein: dict | None = None) -> dict:
    """Bisulfite/EM conversion QC from spike-in controls.
      * Lambda phage DNA is UNmethylated -> residual methylation at any lambda
        cytosine = non-conversion; conversion = 1 - that (all contexts).
      * CpG-methylated pUC19 is fully methylated at CpGs -> 1 - observed CpG
        methylation = over-conversion (false-negative rate)."""
    sp = {**DEFAULT_SPIKEIN, **(spikein or {})}
    cl = meth.chrom.astype(str).str.lower()

    def _mask(names):
        return cl.apply(lambda c: any(c == n or c.startswith(n) for n in names))

    out = {}
    lam = meth[_mask(sp["lambda_contigs"])]
    if len(lam):
        cov = float(np.nansum(lam["cov"]))
        nonconv = float(np.nansum(lam["meth"]) / max(cov, 1))
        out.update(lambda_sites=int(len(lam)),
                   lambda_conversion=round(1 - nonconv, 5),
                   lambda_nonconversion=round(nonconv, 5))
    puc = meth[_mask(sp["puc19_contigs"])]
    pcpg = puc[puc.context == "CpG"] if len(puc) else puc
    if len(pcpg):
        cov = float(np.nansum(pcpg["cov"]))
        methc = float(np.nansum(pcpg["meth"]) / max(cov, 1))
        out.update(puc19_sites=int(len(pcpg)),
                   puc19_methylation=round(methc, 5),
                   puc19_overconversion=round(1 - methc, 5))
    return out


def methylation_qc(meth: pd.DataFrame, gates: dict | None = None,
                   lambda_meth: float | None = None, spikein: dict | None = None) -> dict:
    g = {**DEFAULT_GATES, **(gates or {})}
    # exclude spike-in contigs from the biological CpG metrics
    spk = spikein_conversion(meth, spikein)
    sp = {**DEFAULT_SPIKEIN, **(spikein or {})}
    cl = meth.chrom.astype(str).str.lower()
    is_spike = cl.apply(lambda c: any(c == n or c.startswith(n)
                                      for n in sp["lambda_contigs"] + sp["puc19_contigs"]))
    bio = meth[~is_spike]
    cpg = bio[bio.context == "CpG"]
    chh = bio[bio.context == "CHH"]
    covered = cpg[cpg['cov'] > 0]
    callable_sites = int((cpg['cov'] >= g["min_site_cov"]).sum())

    if "lambda_conversion" in spk:
        conversion = spk["lambda_conversion"]
        conv_source = f"lambda spike-in ({spk['lambda_sites']} cytosines)"
    elif lambda_meth is not None:
        conversion = 1.0 - float(lambda_meth)
        conv_source = "spike-in control (provided value)"
    elif len(chh) > 0:
        conversion = 1.0 - float(np.nansum(chh.meth) / max(np.nansum(chh['cov']), 1))
        conv_source = "non-CpG (CHH) methylation"
    else:
        conversion = np.nan
        conv_source = "not estimable (no spike-in / CHH)"

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
        "spikein": spk,
    }
    metrics["inferred_sex"] = _infer_sex(cpg)

    checks = {
        "conversion": (metrics["conversion_efficiency"] is None) or
                      (metrics["conversion_efficiency"] >= g["min_conversion"]),
        "coverage": metrics["mean_cpg_coverage"] >= g["min_mean_cov"],
        "callable": metrics["n_cpg_callable"] >= g["min_cpg_callable"],
    }
    if "puc19_methylation" in spk:   # over-conversion gate (only when the control is present)
        checks["overconversion"] = spk["puc19_methylation"] >= g["min_puc19_meth"]
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
