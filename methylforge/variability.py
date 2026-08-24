"""Stage 6c — differential variability (DVP / iEVORA).

Complex-disease and autoimmune signals in blood methylation frequently appear as
differences in VARIABILITY between groups rather than in mean methylation
(Teschendorff et al. iEVORA; Paul et al. 2016, T1D discordant twins). This stage
complements the mean-shift DMR test with a variability test.

Per CpG: Bartlett's test for a difference in variance between groups (the
variability signal), combined with a t-test on the means. Following iEVORA, a
differentially variable position (DVP) is called when the variance difference is
significant after BH-FDR AND the means differ at an unadjusted threshold, ranked
by the variance statistic. Needs >= 2 samples per group (>= 3 recommended); with
very small n this is exploratory and hypothesis-generating.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats


def _bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def call_dvps(matrix: pd.DataFrame, groups: dict, *, fdr_var: float = 0.05,
              p_mean: float = 0.05, min_delta_var: float = 0.0) -> pd.DataFrame:
    """iEVORA-style differentially variable positions. matrix: feature x sample beta."""
    cases = [s for s, g in groups.items() if g == "case" and s in matrix.columns]
    ctrls = [s for s, g in groups.items() if g == "control" and s in matrix.columns]
    if len(cases) < 2 or len(ctrls) < 2:
        raise ValueError("DVP analysis needs >= 2 samples per group")
    C = matrix[cases].to_numpy(float)
    K = matrix[ctrls].to_numpy(float)
    n = len(matrix)
    p_var = np.ones(n)
    p_mu = np.ones(n)
    var_c = np.full(n, np.nan)
    var_k = np.full(n, np.nan)
    for i in range(n):
        a = C[i][~np.isnan(C[i])]
        b = K[i][~np.isnan(K[i])]
        if len(a) >= 2 and len(b) >= 2:
            var_c[i] = float(np.var(a, ddof=1))
            var_k[i] = float(np.var(b, ddof=1))
            try:
                if var_c[i] > 0 or var_k[i] > 0:
                    p_var[i] = stats.bartlett(a, b).pvalue
                p_mu[i] = stats.ttest_ind(a, b, equal_var=False).pvalue
            except Exception:
                pass
    q_var = _bh_fdr(p_var)
    out = pd.DataFrame({
        "site": matrix.index,
        "var_case": np.round(var_c, 5), "var_control": np.round(var_k, 5),
        "log2_var_ratio": np.round(np.log2((var_c + 1e-6) / (var_k + 1e-6)), 4),
        "p_var": p_var, "q_var": np.round(q_var, 5), "p_mean": p_mu})
    out["chrom"] = out.site.astype(str).str.split(":").str[0]
    out["pos"] = out.site.astype(str).str.split(":").str[1]
    out["dvp"] = (out.q_var <= fdr_var) & (out.p_mean <= p_mean) & \
                 (np.abs(out.log2_var_ratio) >= min_delta_var)
    out["direction"] = np.where(var_c >= var_k, "hyper-variable(case)", "hypo-variable(case)")
    return out.sort_values("p_var").reset_index(drop=True)


def summarise_dvp(dvps: pd.DataFrame) -> dict:
    n = int(dvps.dvp.sum()) if dvps is not None and "dvp" in dvps else 0
    hv = int(((dvps.dvp) & (dvps.direction.str.startswith("hyper"))).sum()) if n else 0
    return dict(n_dvp=n, n_hypervariable_case=hv,
                headline=(f"{n} differentially variable positions "
                          f"({hv} hyper-variable in cases)." if n
                          else "No differentially variable positions detected."))
