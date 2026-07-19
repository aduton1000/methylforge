"""Stage 8 — region annotation.

Lightweight, dependency-free annotation of DMRs/regions: assigns the nearest
gene from an optional gene BED (chrom, start, end, name), and a genomic-context
label. In the production Nextflow pipeline this stage can instead call `annotatr`
/ a full GTF for CpG-island, promoter, and regulatory-feature context; this
Python implementation keeps the engine runnable without those resources.
"""
from __future__ import annotations
import pandas as pd
import numpy as np


def load_genes_bed(path) -> pd.DataFrame:
    g = pd.read_csv(path, sep="\t", header=None,
                    names=["chrom", "start", "end", "name"], usecols=[0, 1, 2, 3])
    g["chrom"] = g["chrom"].astype(str)
    return g


def annotate_regions(dmrs: pd.DataFrame, genes: pd.DataFrame | None = None) -> pd.DataFrame:
    if dmrs is None or len(dmrs) == 0:
        return dmrs
    out = dmrs.copy()
    if genes is None:
        out["nearest_gene"] = ""
        out["distance"] = np.nan
        return out
    names, dists = [], []
    for _, r in out.iterrows():
        gc = genes[genes.chrom == r.chrom]
        if len(gc) == 0:
            names.append(""); dists.append(np.nan); continue
        mid = (r.start + r.end) / 2
        # distance 0 if overlapping, else gap to nearest edge
        overlap = (gc.start <= r.end) & (gc.end >= r.start)
        if overlap.any():
            names.append(gc[overlap].iloc[0]["name"]); dists.append(0); continue
        d = np.minimum(np.abs(gc.start - mid), np.abs(gc.end - mid))
        i = int(d.values.argmin())
        names.append(gc.iloc[i]["name"]); dists.append(int(d.values[i]))
    out["nearest_gene"] = names
    out["distance"] = dists
    return out
