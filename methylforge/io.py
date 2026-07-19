"""Input/output for MethylForge.

A single normalised in-memory schema for per-cytosine methylation calls, with
readers for the common formats produced by the upstream tools so that every
downstream stage is caller-agnostic:

  normalised columns:  chrom, start, end, context, meth, cov, beta
     - context : 'CpG' | 'CHG' | 'CHH'  (CpG unless stated)
     - meth    : methylated read count (int)
     - cov     : total read count (int)
     - beta    : methylated fraction 0..1  (meth/cov)

Supported inputs (auto-detected by extension / header):
  * Bismark coverage      (.cov / .cov.gz)          chrom start end pct meth unmeth
  * bedGraph methylation  (.bedGraph)               chrom start end pct
  * modkit / ONT bedMethyl(.bed / .bedmethyl)       ENCODE bedMethyl (chrom start end mod score strand ... Nvalid Nmod ...)
  * generic TSV           (has chrom/start/beta cols)
  * methylation-array beta matrix (.beta.tsv)       feature x sample matrix (handled by episignature stage)
"""
from __future__ import annotations
import gzip
import os
import pandas as pd
import numpy as np

REQ = ["chrom", "start", "end", "context", "meth", "cov", "beta"]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _finalise(df: pd.DataFrame) -> pd.DataFrame:
    if "context" not in df:
        df["context"] = "CpG"
    df["cov"] = df["cov"].astype(int)
    df["meth"] = df["meth"].astype(int)
    df["beta"] = np.where(df["cov"] > 0, df["meth"] / df["cov"], np.nan)
    df["chrom"] = df["chrom"].astype(str)
    return df[REQ].sort_values(["chrom", "start"]).reset_index(drop=True)


def read_bismark_cov(path) -> pd.DataFrame:
    df = pd.read_csv(_open(path), sep="\t", header=None,
                     names=["chrom", "start", "end", "pct", "meth", "unmeth"])
    df["cov"] = df["meth"] + df["unmeth"]
    return _finalise(df)


def read_bedmethyl(path) -> pd.DataFrame:
    """ENCODE bedMethyl (modkit / ONT). Columns 10 (Nvalid_cov) and 11 (frac) or Nmod."""
    df = pd.read_csv(_open(path), sep="\t", header=None)
    out = pd.DataFrame({"chrom": df[0], "start": df[1], "end": df[2]})
    # modkit: col9=Nvalid_cov, col10=frac_modified(%), col11=Nmod
    out["cov"] = df[9].astype(int)
    if df.shape[1] > 11:
        out["meth"] = df[11].astype(int)
    else:
        out["meth"] = (df[10].astype(float) / 100.0 * out["cov"]).round().astype(int)
    out["context"] = "CpG"
    return _finalise(out)


def read_bedgraph(path) -> pd.DataFrame:
    df = pd.read_csv(_open(path), sep="\t", header=None,
                     names=["chrom", "start", "end", "pct"], comment="t")
    df["cov"] = 30  # unknown -> nominal; bedGraph carries no coverage
    df["meth"] = (df["pct"] / 100.0 * df["cov"]).round().astype(int)
    return _finalise(df)


def read_methylforge_tsv(path) -> pd.DataFrame:
    """MethylForge native / generic TSV with a header containing chrom,start and (beta or meth+cov)."""
    df = pd.read_csv(_open(path), sep="\t")
    cols = {c.lower(): c for c in df.columns}
    ren = {}
    for k in ["chrom", "start", "end", "context", "meth", "cov", "beta"]:
        if k in cols:
            ren[cols[k]] = k
    df = df.rename(columns=ren)
    if "cov" not in df and "beta" in df:
        df["cov"] = 30
    if "meth" not in df and "beta" in df:
        df["meth"] = (df["beta"] * df["cov"]).round().astype(int)
    if "end" not in df:
        df["end"] = df["start"] + 1
    return _finalise(df)


def read_methylation(path) -> pd.DataFrame:
    """Auto-detect and read any supported per-cytosine methylation file."""
    p = str(path).lower()
    if p.endswith((".cov", ".cov.gz")):
        return read_bismark_cov(path)
    if p.endswith((".bedmethyl", ".bedmethyl.gz")) or ".bedmethyl" in p:
        return read_bedmethyl(path)
    if p.endswith((".bedgraph", ".bedgraph.gz")):
        return read_bedgraph(path)
    # sniff header
    with _open(path) as fh:
        first = fh.readline()
    if "beta" in first.lower() or ("chrom" in first.lower() and "start" in first.lower()):
        return read_methylforge_tsv(path)
    # fall back to bedMethyl (headerless numeric bed)
    return read_bedmethyl(path)


def write_methylation(df: pd.DataFrame, path) -> None:
    df[REQ].to_csv(path, sep="\t", index=False)


def read_dmr_panel(path=None) -> pd.DataFrame:
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "resources", "imprinted_dmrs.tsv")
    df = pd.read_csv(path, sep="\t")
    df["chrom"] = df["chrom"].astype(str)
    return df
