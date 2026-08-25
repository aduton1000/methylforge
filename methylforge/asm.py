"""Stage 6b — allele-specific methylation (ASM).

Detects methylation that differs between the two alleles/haplotypes — the basis
of genomic imprinting and of cis-regulatory allelic effects (e.g. a haplotype
carrying a disease variant regulated differently from its wild-type counterpart).
Three entry points:

  * Haplotype-resolved (rigorous): two per-cytosine methylation files, one per
    haplotype — e.g. `modkit pileup --partition-tag HP` (long reads) or a
    bisulfite BAM split by a phased heterozygous SNP. Per CpG, a Fisher exact
    test on (methylated, unmethylated) counts between haplotypes, BH-FDR, and an
    effect-size (Δβ) filter -> ASM sites -> merged ASM regions.
  * Long-format allele table: chrom,pos,allele,meth,cov with >= 2 alleles ->
    pivoted to the two best-covered alleles, then as above.
  * Single-sample candidate scan (no phasing): flags runs of CpGs held at
    intermediate methylation (~0.5) with adequate coverage as *candidate* ASM
    (imprinting-like). A heuristic for prioritising loci (e.g. around a gene of
    interest) that MUST be confirmed with phased/haplotype data — aggregate short
    reads cannot by themselves distinguish true allelic 0/1 splits from uniform
    intermediate methylation.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from . import io as eio


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


def _prep(meth: pd.DataFrame, min_cov: int) -> pd.DataFrame:
    m = meth[(meth.context == "CpG") & (meth["cov"] >= min_cov)].copy()
    m["unmeth"] = (m["cov"] - m["meth"]).clip(lower=0)
    m["key"] = m.chrom.astype(str) + ":" + m.start.astype(int).astype(str)
    return m[["key", "chrom", "start", "meth", "unmeth", "cov", "beta"]]


def asm_from_haplotypes(hp1: pd.DataFrame, hp2: pd.DataFrame, *, min_cov: int = 5,
                        min_delta: float = 0.20, fdr: float = 0.05,
                        hp1_name: str = "hp1", hp2_name: str = "hp2") -> pd.DataFrame:
    """Per-CpG allele-specific methylation between two haplotype methylation frames."""
    a = _prep(hp1, min_cov).set_index("key")
    b = _prep(hp2, min_cov).set_index("key")
    shared = a.index.intersection(b.index)
    if len(shared) == 0:
        return pd.DataFrame(columns=["site", "chrom", "pos", f"beta_{hp1_name}", f"beta_{hp2_name}",
                                     "delta", "pvalue", "qvalue", "asm", "direction"])
    a = a.loc[shared]
    b = b.loc[shared]
    pvals = np.ones(len(shared))
    for i, k in enumerate(shared):
        table = [[int(a.at[k, "meth"]), int(a.at[k, "unmeth"])],
                 [int(b.at[k, "meth"]), int(b.at[k, "unmeth"])]]
        try:
            pvals[i] = stats.fisher_exact(table)[1]
        except Exception:
            pvals[i] = 1.0
    delta = a["beta"].to_numpy(float) - b["beta"].to_numpy(float)
    q = _bh_fdr(pvals)
    out = pd.DataFrame({
        "site": list(shared),
        "chrom": a["chrom"].to_numpy(),
        "pos": a["start"].to_numpy(int),
        f"beta_{hp1_name}": np.round(a["beta"].to_numpy(float), 4),
        f"beta_{hp2_name}": np.round(b["beta"].to_numpy(float), 4),
        "delta": np.round(delta, 4),
        "pvalue": pvals, "qvalue": np.round(q, 5)})
    out["asm"] = (np.abs(out.delta) >= min_delta) & (out.qvalue <= fdr)
    out["direction"] = np.where(out.delta > 0, f"{hp1_name}_hyper", f"{hp2_name}_hyper")
    return out.sort_values(["chrom", "pos"]).reset_index(drop=True)


def merge_asm_regions(sites: pd.DataFrame, max_gap: int = 1000, min_cpg: int = 3) -> pd.DataFrame:
    """Merge adjacent, same-direction ASM CpGs into ASM regions."""
    if sites is None or len(sites) == 0 or "asm" not in sites:
        return pd.DataFrame(columns=["chrom", "start", "end", "n_cpg", "direction",
                                     "mean_delta", "mean_abs_delta", "min_qvalue"])
    sig = sites[sites.asm].copy()
    regions = []
    for (chrom, direction), grp in sig.groupby(["chrom", "direction"]):
        grp = grp.sort_values("pos")
        buf, prev = [], None
        for _, row in grp.iterrows():
            if prev is None or row.pos - prev <= max_gap:
                buf.append(row)
            else:
                regions.append(_finalise(chrom, direction, buf))
                buf = [row]
            prev = row.pos
        if buf:
            regions.append(_finalise(chrom, direction, buf))
    reg = pd.DataFrame(regions)
    if len(reg):
        reg = reg[reg.n_cpg >= min_cpg].sort_values("mean_abs_delta", ascending=False).reset_index(drop=True)
    return reg


def _finalise(chrom, direction, buf):
    rows = pd.DataFrame(buf)
    return dict(chrom=chrom, start=int(rows.pos.min()), end=int(rows.pos.max()),
                width=int(rows.pos.max() - rows.pos.min() + 1),
                n_cpg=int(len(rows)), direction=direction,
                mean_delta=round(float(rows.delta.mean()), 4),
                mean_abs_delta=round(float(rows.delta.abs().mean()), 4),
                min_qvalue=float(rows.qvalue.min()))


def read_allele_table(path) -> tuple[pd.DataFrame, pd.DataFrame, tuple]:
    """Long-format allele table (chrom,pos,allele,meth,cov) -> the two best-covered
    alleles as (hp1_meth, hp2_meth, (name1, name2)) in the normalised io schema."""
    df = pd.read_csv(path, sep=None, engine="python")
    lc = {c.lower(): c for c in df.columns}
    ch, po, al = lc["chrom"], lc.get("pos", lc.get("start")), lc["allele"]
    me, co = lc["meth"], lc["cov"]
    cover = df.groupby(df[al])[co].sum().sort_values(ascending=False)
    if len(cover) < 2:
        raise ValueError("allele table needs >= 2 distinct alleles")
    a1, a2 = cover.index[0], cover.index[1]

    def _frame(allele):
        s = df[df[al] == allele]
        out = pd.DataFrame({"chrom": s[ch].astype(str), "start": s[po].astype(int)})
        out["end"] = out["start"] + 1
        out["context"] = "CpG"
        out["meth"] = s[me].astype(int).values
        out["cov"] = s[co].astype(int).values
        out["beta"] = np.where(out["cov"] > 0, out["meth"] / out["cov"], np.nan)
        return out.reset_index(drop=True)

    return _frame(a1), _frame(a2), (str(a1), str(a2))


def candidate_asm_scan(meth: pd.DataFrame, *, min_cov: int = 10, lo: float = 0.35,
                       hi: float = 0.65, min_sites: int = 4, max_gap: int = 1000) -> pd.DataFrame:
    """HEURISTIC single-sample scan (no phasing): runs of adequately-covered CpGs
    held at intermediate methylation (lo..hi) — candidate ASM / imprinting-like
    regions to prioritise for phased confirmation. Not a statistical ASM call."""
    m = meth[(meth.context == "CpG") & (meth["cov"] >= min_cov) &
             (meth.beta >= lo) & (meth.beta <= hi)].copy()
    if len(m) == 0:
        return pd.DataFrame(columns=["chrom", "start", "end", "n_cpg", "mean_beta", "min_cov"])
    m = m.sort_values(["chrom", "start"])
    regions, buf, prev, cur_chrom = [], [], None, None
    for _, row in m.iterrows():
        if cur_chrom != row.chrom or (prev is not None and row.start - prev > max_gap):
            if len(buf) >= min_sites:
                regions.append(_finalise_candidate(buf))
            buf = []
        buf.append(row)
        prev, cur_chrom = row.start, row.chrom
    if len(buf) >= min_sites:
        regions.append(_finalise_candidate(buf))
    return pd.DataFrame(regions)


def _finalise_candidate(buf):
    rows = pd.DataFrame(buf)
    return dict(chrom=str(rows.chrom.iloc[0]), start=int(rows.start.min()),
                end=int(rows.start.max()), width=int(rows.start.max() - rows.start.min() + 1),
                n_cpg=int(len(rows)), mean_beta=round(float(rows.beta.mean()), 4),
                min_cov=int(rows["cov"].min()), note="candidate ASM (heuristic; confirm with phasing)")


def summarise_asm(sites: pd.DataFrame, regions: pd.DataFrame) -> dict:
    n_sites = int(sites.asm.sum()) if sites is not None and "asm" in sites else 0
    n_reg = int(len(regions)) if regions is not None else 0
    return dict(n_asm_cpg=n_sites, n_asm_regions=n_reg,
                headline=(f"{n_reg} allele-specific methylation region(s) from {n_sites} ASM CpGs."
                          if n_reg else "No significant allele-specific methylation detected."))


def load_pair(hp1_path, hp2_path):
    """Convenience: read two haplotype methylation files in any supported format."""
    return eio.read_methylation(hp1_path), eio.read_methylation(hp2_path)
