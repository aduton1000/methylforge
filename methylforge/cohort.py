"""Cohort matrix builder.

Assembles a feature x sample methylation (beta) matrix from many per-sample
methylation files (any format `io.read_methylation` accepts). The matrix is the
input to the DMR stage (case-vs-control) and to episignature reference building /
classification.

Design choices (best-practice):
  * per-site coverage filter (min_cov) applied to each sample before merging;
  * a site is kept only if covered in >= min_samples_frac of samples (avoids
    columns of mostly-missing data);
  * optional restriction to a feature list (e.g. episignature CpGs) or a region
    BED (targeted panel);
  * context-aware (CpG by default).
"""
from __future__ import annotations
import os
import math
import glob
import numpy as np
import pandas as pd
from . import io as eio


def _load_sites(path) -> list[str]:
    df = pd.read_csv(path, sep=None, engine="python", header=None)
    if df.shape[1] == 1:
        return df.iloc[:, 0].astype(str).tolist()
    return [f"{a}:{int(b)}" for a, b in zip(df.iloc[:, 0], df.iloc[:, 1])]


def _load_bed(path) -> pd.DataFrame:
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                    names=["chrom", "start", "end"])
    b["chrom"] = b["chrom"].astype(str)
    return b


def _in_regions(chroms: np.ndarray, starts: np.ndarray, bed: pd.DataFrame) -> np.ndarray:
    keep = np.zeros(len(chroms), bool)
    for chrom, grp in bed.groupby("chrom"):
        idx = np.where(chroms == chrom)[0]
        if len(idx) == 0:
            continue
        s = starts[idx]
        m = np.zeros(len(idx), bool)
        for a, e in zip(grp.start.values, grp.end.values):
            m |= (s >= a) & (s < e)
        keep[idx] = m
    return keep


def build_matrix(sample_paths: dict, min_cov: int = 5, min_samples_frac: float = 0.75,
                 sites=None, region_bed: str | None = None, context: str = "CpG") -> pd.DataFrame:
    bed = _load_bed(region_bed) if region_bed else None
    per = {}
    for name, path in sample_paths.items():
        m = eio.read_methylation(path)
        m = m[(m.context == context) & (m["cov"] >= min_cov)]
        if bed is not None:
            m = m[_in_regions(m.chrom.values.astype(str), m.start.values, bed)]
        ids = [f"{c}:{p}" for c, p in zip(m.chrom, m.start)]
        s = pd.Series(m.beta.values, index=ids)
        per[name] = s[~s.index.duplicated(keep="first")]
    mat = pd.DataFrame(per)
    if sites is not None:
        mat = mat.reindex(sites)
    n = mat.shape[1]
    min_n = max(1, math.ceil(min_samples_frac * n))
    mat = mat[mat.notna().sum(axis=1) >= min_n]
    mat.index.name = "feature"
    return mat.round(4)


def read_samplesheet(path: str) -> tuple[dict, dict]:
    """CSV with columns: sample, (methylation|path|file), optional (group|class).
    Relative paths are resolved against the samplesheet's directory."""
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    scol = cols.get("sample")
    pcol = cols.get("methylation") or cols.get("path") or cols.get("file")
    gcol = cols.get("group") or cols.get("class")
    if not scol or not pcol:
        raise ValueError("samplesheet needs at least 'sample' and 'methylation' (or path/file) columns")
    base = os.path.dirname(os.path.abspath(path))
    paths, groups = {}, {}
    for _, r in df.iterrows():
        p = str(r[pcol])
        if not os.path.isabs(p):
            cand = os.path.join(base, p)
            p = cand if os.path.isfile(cand) else p
        paths[str(r[scol])] = p
        if gcol and pd.notna(r[gcol]):
            groups[str(r[scol])] = str(r[gcol])
    return paths, groups


def paths_from_globs(patterns) -> dict:
    paths = {}
    for pat in patterns:
        for f in sorted(glob.glob(pat)):
            name = os.path.basename(f).split(".")[0]
            paths[name] = f
    return paths


def write_groups(groups: dict, path: str) -> None:
    pd.DataFrame({"sample": list(groups), "group": list(groups.values())}).to_csv(
        path, sep="\t", index=False)
