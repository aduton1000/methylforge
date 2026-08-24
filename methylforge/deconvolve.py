"""Stage 4b — reference-based cell-type deconvolution.

Estimates the cell-type composition of a bulk methylation profile (e.g. whole
blood) from a reference of cell-type-specific methylation signatures, so that
composition can be (i) reported as a result in its own right and (ii) carried as
a covariate into the DMR model. This is essential whenever cases and controls
differ in cell make-up — the classic confounder in blood methylation studies
(e.g. paediatric vs adult blood, or an immune-skewed patient vs a healthy
relative) — because uncorrected composition differences masquerade as disease
DMRs.

Two estimators (both pure numpy/scipy, no sklearn dependency):
  * 'CP'  — constrained projection (Houseman et al. 2012): non-negative
            proportions that sum to 1, least-squares fit. Default; robust for
            small marker panels.
  * 'RPC' — robust partial correlations (Teschendorff/EpiDISH-style): Huber IRLS
            regression, coefficients clipped at 0 and renormalised; less
            sensitive to a few outlying markers.

Reference format (TSV): first column = feature id ('chrom:pos' for sequencing,
or a probe id for arrays), one column per cell type = that cell type's mean
methylation beta (0..1). Supply a *validated* reference — e.g. the Salas et al.
2022 12-cell IDOL panel, coordinate-mapped to your genome for WGBS/EM-seq, or the
IDOL 450K/EPIC panel for array data. Format + sourcing notes are in
resources/README_references.md. No reference is redistributed with MethylForge.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.optimize import nnls


def load_reference(path) -> pd.DataFrame:
    """Load a feature x cell-type reference (beta 0..1)."""
    ref = pd.read_csv(path, sep="\t", index_col=0)
    ref.index = ref.index.astype(str)
    ref.index.name = "feature"
    return ref


def _align(query: pd.Series, ref: pd.DataFrame):
    q = query.copy()
    q.index = q.index.astype(str)
    q = q[~q.index.duplicated(keep="first")]
    shared = [f for f in ref.index if f in q.index and np.isfinite(q.get(f, np.nan))]
    need = max(2, ref.shape[1])
    if len(shared) < need:
        raise ValueError(f"only {len(shared)} shared finite markers between query and reference "
                         f"(need >= {need}); check the feature-id keyspace (chrom:pos vs probe id)")
    R = ref.loc[shared].to_numpy(float)
    y = q.loc[shared].to_numpy(float)
    return R, y, shared


def _cp(R, y):
    """Non-negative, sum-to-one least squares (Houseman constrained projection):
    NNLS gives non-negativity; sum-to-one is imposed by renormalising."""
    w, _ = nnls(R, y)
    s = w.sum()
    return w / s if s > 0 else w


def _rpc(R, y, n_iter=50, c=1.345, tol=1e-6):
    """Huber IRLS robust regression; coefficients clipped >= 0 and renormalised."""
    n = R.shape[0]
    w = np.ones(n)
    beta = np.zeros(R.shape[1])
    for _ in range(n_iter):
        sw = np.sqrt(w)
        beta_new, *_ = np.linalg.lstsq(R * sw[:, None], y * sw, rcond=None)
        resid = y - R @ beta_new
        mad = np.median(np.abs(resid - np.median(resid))) * 1.4826
        s = mad if mad > 0 else (np.std(resid) or 1.0)
        r = np.abs(resid) / s
        w = np.where(r <= c, 1.0, c / np.maximum(r, 1e-9))
        step = np.max(np.abs(beta_new - beta))
        beta = beta_new
        if step < tol:
            break
    beta = np.clip(beta, 0, None)
    s = beta.sum()
    return beta / s if s > 0 else beta


def deconvolve(query: pd.Series, ref: pd.DataFrame, method: str = "CP") -> dict:
    """Estimate cell-type proportions for one sample. query: feature -> beta."""
    R, y, shared = _align(query, ref)
    w = _rpc(R, y) if method.upper() == "RPC" else _cp(R, y)
    fitted = R @ w
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return dict(proportions={ct: round(float(v), 4) for ct, v in zip(ref.columns, w)},
                method=method.upper(), n_markers=len(shared),
                rmse=round(float(np.sqrt(np.mean((y - fitted) ** 2))), 4),
                fit_r2=round(float(r2), 4) if r2 == r2 else None)


def deconvolve_matrix(matrix: pd.DataFrame, ref: pd.DataFrame, method: str = "CP") -> pd.DataFrame:
    """Deconvolve every sample (column) of a feature x sample beta matrix ->
    a samples x cell-types proportion table (NaN row if a sample lacks markers)."""
    rows = {}
    for s in matrix.columns:
        try:
            rows[s] = deconvolve(matrix[s], ref, method=method)["proportions"]
        except ValueError:
            rows[s] = {ct: np.nan for ct in ref.columns}
    return pd.DataFrame(rows).T.reindex(columns=list(ref.columns))


def summarise(result: dict) -> dict:
    props = result["proportions"]
    top = max(props.items(), key=lambda kv: kv[1]) if props else ("", 0.0)
    return dict(headline=f"Estimated composition ({result['method']}, {result['n_markers']} markers): "
                         + ", ".join(f"{k} {v:.1%}" for k, v in
                                     sorted(props.items(), key=lambda kv: -kv[1])),
                dominant=top[0], dominant_fraction=round(float(top[1]), 4),
                fit_r2=result.get("fit_r2"))
