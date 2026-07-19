"""Illumina Infinium MethylationEPIC (850K) / 450K array ingestion.

The IDAT -> beta-matrix step (minfi/sesame QC + normalisation) runs in R via the
Nextflow module `modules/stage0_array.nf` / `bin/epic_process.R` (the R/Bioconductor
env `env/array.yml`). This module takes the resulting **probe x sample beta matrix**
and maps probes to genomic coordinates using an EPIC manifest, producing the SAME
coordinate-keyed schema (`chrom:pos`) used everywhere else in MethylForge — so array
data feeds the imprinting, cohort, episignature and DMR stages unchanged.

A user supplies the manifest (from minfi's IlluminaHumanMethylationEPICanno, a TSV
with probe/Name, chr, pos columns) — MethylForge does not bundle the licensed manifest.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd


def read_beta_matrix(path) -> pd.DataFrame:
    """Probe/feature x sample beta matrix; first column is the probe/feature id."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df


def _col(df, *names):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def load_manifest(path) -> pd.DataFrame:
    """EPIC manifest TSV -> probe, chrom, pos. Tolerant of column naming
    (Name/IlmnID/probe; chr/chrom/CHR; pos/MAPINFO/position)."""
    df = pd.read_csv(path, sep="\t")
    pcol = _col(df, "probe", "Name", "IlmnID", "probeID")
    ccol = _col(df, "chrom", "chr", "CHR")
    poscol = _col(df, "pos", "MAPINFO", "position", "start")
    if not (pcol and ccol and poscol):
        raise ValueError("manifest needs probe (Name/IlmnID), chr, and pos columns")
    out = pd.DataFrame({"probe": df[pcol].astype(str),
                        "chrom": df[ccol].astype(str),
                        "pos": df[poscol].astype("Int64")}).dropna()
    out["chrom"] = out["chrom"].apply(lambda c: c if str(c).startswith("chr") else f"chr{c}")
    return out


def map_probes_to_coords(beta: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """Map a probe x sample beta matrix to a coord (chrom:pos) x sample matrix."""
    m = manifest.set_index("probe")
    common = beta.index.intersection(m.index)
    b = beta.loc[common].copy()
    ids = [f"{m.loc[p, 'chrom']}:{int(m.loc[p, 'pos'])}" for p in common]
    b.index = ids
    b = b.groupby(level=0).mean()          # collapse any duplicate positions
    b.index.name = "feature"
    return b.round(4)


def to_sample_methyl(coord_matrix: pd.DataFrame, sample: str, out_path: str,
                     cov: int = 30, context: str = "CpG") -> str:
    """Write one array sample as a MethylForge per-cytosine methylation TSV so the
    single-sample `run` (QC + imprinting + episignature) works on array data.
    Arrays carry no coverage, so a nominal `cov` is assigned (documented in QC)."""
    s = coord_matrix[sample].dropna()
    chrom = [i.split(":")[0] for i in s.index]
    pos = [int(i.split(":")[1]) for i in s.index]
    meth = np.round(s.values * cov).astype(int)
    df = pd.DataFrame({"chrom": chrom, "start": pos, "end": [p + 1 for p in pos],
                       "context": context, "meth": meth, "cov": cov,
                       "beta": s.values.round(4)})
    df.sort_values(["chrom", "start"]).to_csv(out_path, sep="\t", index=False)
    return out_path


def ingest(beta_path, manifest_path=None, out_matrix=None, out_methyl_dir=None,
           cov=30) -> dict:
    beta = read_beta_matrix(beta_path)
    if manifest_path:
        coord = map_probes_to_coords(beta, load_manifest(manifest_path))
    else:
        coord = beta  # already coordinate-keyed
        coord.index.name = "feature"
    result = {"n_features": coord.shape[0], "n_samples": coord.shape[1]}
    if out_matrix:
        coord.to_csv(out_matrix, sep="\t")
        result["matrix"] = out_matrix
    if out_methyl_dir:
        os.makedirs(out_methyl_dir, exist_ok=True)
        files = {}
        for s in coord.columns:
            p = os.path.join(out_methyl_dir, f"{s}.methyl.tsv")
            to_sample_methyl(coord, s, p, cov=cov)
            files[s] = p
        result["methyl_files"] = files
    return result
