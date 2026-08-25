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


def _mvalue(x):
    """logit / M-value transform (respects the bounded 0..1 beta scale; the standard
    for parametric methylation differential testing, cf. limma/methylKit)."""
    x = np.clip(x, 1e-3, 1 - 1e-3)
    return np.log2(x / (1 - x))


def call_dmps(matrix: pd.DataFrame, groups: dict, min_delta=0.10,
              fdr=0.05, test="auto") -> pd.DataFrame:
    """test: 'auto' (Welch t-test on M-values for cohorts <12/group, where the
    Mann-Whitney asymptotic p-floor is too high to survive FDR; the Mann-Whitney
    rank test for larger cohorts), 'ttest', or 'mannwhitney'."""
    cases = [s for s, g in groups.items() if g == "case" and s in matrix.columns]
    ctrls = [s for s, g in groups.items() if g == "control" and s in matrix.columns]
    if not cases or not ctrls:
        raise ValueError("need >=1 case and >=1 control sample present in the matrix")
    C = matrix[cases].to_numpy(float)
    K = matrix[ctrls].to_numpy(float)
    delta = np.nanmean(C, axis=1) - np.nanmean(K, axis=1)
    use = test
    if test == "auto":
        use = "mannwhitney" if min(len(cases), len(ctrls)) >= 12 else "ttest"
    pvals = np.ones(len(matrix))
    for i in range(len(matrix)):
        a, b = C[i][~np.isnan(C[i])], K[i][~np.isnan(K[i])]
        if len(a) >= 1 and len(b) >= 1 and (len(a) + len(b)) >= 3:
            try:
                if use == "mannwhitney" and len(a) >= 2 and len(b) >= 2:
                    pvals[i] = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
                else:
                    pvals[i] = stats.ttest_ind(_mvalue(a), _mvalue(b), equal_var=False).pvalue
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


SEX_CHROMS = {"X", "Y", "M", "MT", "chrX", "chrY", "chrM", "chrMT"}


def restrict_autosomes(matrix: pd.DataFrame) -> pd.DataFrame:
    """Drop chrX/chrY/chrM rows from a feature x sample matrix (index 'chrom:pos').
    Essential before an affected-vs-unaffected contrast on mixed-sex samples, where
    X-inactivation otherwise dominates the result."""
    chrom = matrix.index.astype(str).str.split(":").str[0].str.replace("chr", "", regex=False)
    return matrix[~chrom.isin({"X", "Y", "M", "MT"})]


def _build_design(samples, groups, covariates=None, pair=None):
    """Design matrix for the covariate/paired DMR model: intercept + group(case=1)
    + numeric/one-hot covariates + one-hot pair(block) effects. Returns
    (X DataFrame indexed by sample, name of the group column)."""
    X = pd.DataFrame(index=list(samples))
    X["(Intercept)"] = 1.0
    X["group"] = [1.0 if groups.get(s) == "case" else 0.0 for s in samples]
    if covariates is not None:
        cov = covariates.reindex(samples)
        for c in cov.columns:
            col = cov[c]
            if pd.api.types.is_numeric_dtype(col):
                X[c] = col.astype(float).to_numpy()
            else:
                d = pd.get_dummies(col.astype(str), prefix=c, drop_first=True)
                for dc in d.columns:
                    X[dc] = d[dc].astype(float).to_numpy()
    if pair is not None:
        ps = pd.Series({s: str(pair.get(s, "NA")) for s in samples}).reindex(samples)
        d = pd.get_dummies(ps, prefix="pair", drop_first=True)
        for dc in d.columns:
            X[dc] = d[dc].astype(float).to_numpy()
    return X, "group"


def call_dmps_glm(matrix: pd.DataFrame, groups: dict, covariates: pd.DataFrame | None = None,
                  pair: dict | None = None, min_delta=0.10, fdr=0.05) -> pd.DataFrame:
    """Covariate-adjusted / paired differential methylation. Fits, per CpG, an OLS
    on M-values:  M ~ group + covariates (+ pair block), and tests the group
    coefficient — so sex, age and cell-composition can be adjusted out, and a twin
    pair can be modelled as a block. Output schema matches `call_dmps` (so
    `merge_dmrs` applies). Δβ is reported in beta space for interpretability.

    covariates: a samples x covariates DataFrame (numeric or categorical).
    pair: {sample: pair_id} fixed-effect blocking (e.g. twin/family)."""
    samples = [s for s in matrix.columns if s in groups]
    cases = [s for s in samples if groups[s] == "case"]
    ctrls = [s for s in samples if groups[s] == "control"]
    if not cases or not ctrls:
        raise ValueError("need >=1 case and >=1 control sample present in the matrix")
    X, gcol = _build_design(samples, groups, covariates, pair)
    Xmat = X.to_numpy(float)
    gidx = list(X.columns).index(gcol)
    p = Xmat.shape[1]

    C = matrix[cases].to_numpy(float)
    K = matrix[ctrls].to_numpy(float)
    delta = np.nanmean(C, axis=1) - np.nanmean(K, axis=1)
    Y = _mvalue(matrix[samples].to_numpy(float))   # sites x samples (may contain NaN)
    nsite = Y.shape[0]
    pvals = np.ones(nsite)
    coefs = np.zeros(nsite)

    complete = ~np.isnan(Y).any(axis=1)
    # Fast vectorised path for rows with no missing samples
    if complete.any() and Xmat.shape[0] - p > 0:
        try:
            XtX_inv = np.linalg.inv(Xmat.T @ Xmat)
            H = XtX_inv @ Xmat.T                       # p x n
            Yc = Y[complete]                           # m x n
            B = Yc @ H.T                               # m x p
            resid = Yc - B @ Xmat.T                    # m x n
            dof = Xmat.shape[0] - p
            sigma2 = (resid ** 2).sum(axis=1) / dof
            se = np.sqrt(sigma2 * XtX_inv[gidx, gidx])
            with np.errstate(divide="ignore", invalid="ignore"):
                t = np.where(se > 0, B[:, gidx] / se, 0.0)
            pc = 2 * stats.t.sf(np.abs(t), dof)
            idx = np.where(complete)[0]
            pvals[idx] = np.where(se > 0, pc, 1.0)
            coefs[idx] = B[:, gidx]
        except np.linalg.LinAlgError:
            complete[:] = False   # fall through to per-row loop
    # Per-row path for rows with missing data
    for i in np.where(~complete)[0]:
        yi = Y[i]
        ok = ~np.isnan(yi)
        if ok.sum() - p <= 0:
            continue
        Xi = Xmat[ok]
        try:
            XtX_inv = np.linalg.inv(Xi.T @ Xi)
        except np.linalg.LinAlgError:
            continue
        beta = XtX_inv @ Xi.T @ yi[ok]
        resid = yi[ok] - Xi @ beta
        dof = ok.sum() - p
        sigma2 = float(resid @ resid) / dof
        se = np.sqrt(sigma2 * XtX_inv[gidx, gidx])
        if se > 0:
            tval = beta[gidx] / se
            pvals[i] = 2 * stats.t.sf(abs(tval), dof)
            coefs[i] = beta[gidx]

    q = _bh_fdr(pvals)
    out = pd.DataFrame({
        "site": matrix.index, "delta": np.round(delta, 4),
        "coef_M": np.round(coefs, 4),
        "case_mean": np.round(np.nanmean(C, axis=1), 4),
        "control_mean": np.round(np.nanmean(K, axis=1), 4),
        "pvalue": pvals, "qvalue": np.round(q, 5)})
    out["chrom"] = out.site.astype(str).str.split(":").str[0]
    out["pos"] = out.site.astype(str).str.split(":").str[1].astype(int)
    out["dmp"] = (np.abs(out.delta) >= min_delta) & (out.qvalue <= fdr)
    out["direction"] = np.where(out.delta > 0, "hyper", "hypo")
    return out.sort_values(["chrom", "pos"]).reset_index(drop=True)


def summarise_dmr(dmps: pd.DataFrame, dmrs: pd.DataFrame) -> dict:
    n_dmp = int(dmps.dmp.sum())
    n_dmr = int(len(dmrs))
    hyper = int((dmrs.direction == "hyper").sum()) if n_dmr else 0
    hypo = int((dmrs.direction == "hypo").sum()) if n_dmr else 0
    return dict(n_tested=int(len(dmps)), n_dmp=n_dmp, n_dmr=n_dmr,
                n_hyper=hyper, n_hypo=hypo,
                headline=f"{n_dmr} DMRs ({hyper} hyper-, {hypo} hypomethylated) "
                         f"from {n_dmp} significant CpGs.")
