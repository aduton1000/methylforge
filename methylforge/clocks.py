"""Stage 4c — epigenetic clocks (DNAm age & age acceleration).

A generic linear epigenetic-clock engine: DNAmAge = transform(intercept + Σ coef_i·β_i).
It covers purely additive clocks (Hannum, PhenoAge, DunedinPACE-style pace models)
and the Horvath pan-tissue clock with its calibrated non-linear age transform.
Clocks are per-sample, so they work at any n — including this study's small twin
design.

Coefficient tables are NOT redistributed (they belong to their publications).
Supply each clock as a TSV with a `feature` and a `coefficient` column, one row
per marker plus a special `(Intercept)` row; an optional `mean` column supplies
per-marker gold-standard means used to impute markers missing from the sample. A
sidecar `<clock>.meta.json` (or `meta_path`) sets the transform:
    {"name": "Horvath2013", "transform": "horvath", "adult_age": 20}
Feature ids must match the query keyspace: use 'chrom:pos' for sequencing data
(map array probe ids with `methylforge array`), or probe ids for array betas.

Intrinsic vs extrinsic acceleration: pass cell-type proportions (from the
`deconvolve` stage) to regress DNAmAge on composition and report the
composition-independent (intrinsic) residual — the immune-cell-independent signal.
"""
from __future__ import annotations
import json
import os
import numpy as np
import pandas as pd


def load_clock(path, meta_path: str | None = None) -> dict:
    """Load a clock coefficient table (+ optional meta) into a clock dict."""
    df = pd.read_csv(path, sep="\t")
    lc = {c.lower(): c for c in df.columns}
    fcol = lc.get("feature", df.columns[0])
    ccol = lc.get("coefficient", lc.get("coef", df.columns[1]))
    feats = df[fcol].astype(str)
    is_int = feats.str.lower().isin(["(intercept)", "intercept"])
    intercept = float(df.loc[is_int, ccol].iloc[0]) if is_int.any() else 0.0
    body = df[~is_int]
    coef = pd.Series(body[ccol].to_numpy(float), index=body[fcol].astype(str).to_numpy())
    means = None
    if "mean" in lc:
        means = pd.Series(body[lc["mean"]].to_numpy(float), index=body[fcol].astype(str).to_numpy())
    meta = {"name": os.path.basename(str(path)), "transform": "none", "adult_age": 20.0}
    if meta_path is None:
        cand = str(path) + ".meta.json"
        meta_path = cand if os.path.isfile(cand) else None
    if meta_path and os.path.isfile(meta_path):
        meta.update(json.load(open(meta_path)))
    return dict(intercept=intercept, coef=coef, means=means, meta=meta)


def clock_from_series(coef: pd.Series, intercept: float = 0.0, name: str = "clock",
                      transform: str = "none", adult_age: float = 20.0,
                      means: pd.Series | None = None) -> dict:
    """Build a clock dict programmatically (used by tests / in-memory clocks)."""
    return dict(intercept=float(intercept), coef=coef.astype(float), means=means,
                meta={"name": name, "transform": transform, "adult_age": adult_age})


def _horvath(x, adult=20.0):
    return (1 + adult) * np.exp(x) - 1 if x < 0 else (1 + adult) * x + adult


def predict_age(query: pd.Series, clock: dict) -> dict:
    """Predict DNAm age for one sample. query: feature -> beta (0..1)."""
    q = query.copy()
    q.index = q.index.astype(str)
    q = q[~q.index.duplicated(keep="first")]
    coef = clock["coef"]
    feats = coef.index
    beta = q.reindex(feats)
    n_missing = int(beta.isna().sum())
    if clock.get("means") is not None:
        beta = beta.fillna(clock["means"].reindex(feats))
    beta = beta.fillna(0.5)   # last-resort neutral impute
    raw = float(clock["intercept"] + float(np.dot(coef.to_numpy(float), beta.to_numpy(float))))
    tr = clock["meta"].get("transform", "none")
    age = _horvath(raw, clock["meta"].get("adult_age", 20.0)) if tr == "horvath" else raw
    return dict(clock=clock["meta"].get("name", "clock"), transform=tr,
                dnam_age=round(float(age), 3), linear_predictor=round(raw, 5),
                n_clock_features=int(len(feats)), n_missing=n_missing,
                frac_missing=round(n_missing / max(len(feats), 1), 3))


def run_clocks(query: pd.Series, clocks: dict) -> pd.DataFrame:
    """Apply several clocks to one sample. clocks: {name: clock_dict}."""
    return pd.DataFrame([{**predict_age(query, c)} for c in clocks.values()])


def delta_age(dnam_age, chrono_age) -> float:
    """Simple age acceleration for a single sample (DNAmAge - chronological age)."""
    return round(float(dnam_age) - float(chrono_age), 3)


def cohort_acceleration(dnam_ages, chrono_ages) -> np.ndarray:
    """Residual of DNAmAge ~ chronological age across a cohort (>= 3 samples) —
    the standard 'age acceleration' definition."""
    x = np.asarray(chrono_ages, float)
    y = np.asarray(dnam_ages, float)
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ coef


def intrinsic_acceleration(dnam_ages, proportions: pd.DataFrame) -> np.ndarray:
    """Composition-independent (intrinsic) age acceleration: residual of DNAmAge on
    cell-type proportions across samples. proportions rows align to dnam_ages order."""
    y = np.asarray(dnam_ages, float)
    X = np.column_stack([np.ones(len(y)), np.asarray(proportions, float)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ coef
