"""Stage 2b — sample-identity / swap QC.

First-line detection of sample swaps and unexpected duplicates directly from the
methylation data — the check that catches a mislabelled tube before it becomes a
false biological finding:

  * pairwise correlation of genome-wide beta profiles across a cohort (a true
    duplicate/relabel ~1; a swap makes the "wrong" pair the most similar);
  * a nearest-neighbour table + duplicate/outlier flags;
  * a coverage-based sex estimate per sample (chrY/autosome ratio) to cross-check
    the declared sex — the same estimator used in the QC stage.

Full genotype concordance to an external call set (e.g. confirming these
methylomes are the same individuals as the WES samples) needs SNPs from the
alignment — call genotypes from the bisulfite BAM (allowing C/T ambiguity on the
converted strand) and compare to the exome VCF. That needs the BAM, so it lives
in the Nextflow layer; this module covers the matrix-only checks that need no
extra data.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from . import io as eio, qc as eqc


def correlation_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    """Pairwise Pearson correlation of sample beta profiles (pairwise-complete)."""
    return matrix.corr(method="pearson", min_periods=1)


def nearest_neighbours(corr: pd.DataFrame) -> pd.DataFrame:
    """For each sample, its most-correlated other sample and that correlation."""
    rows = []
    for s in corr.columns:
        others = corr[s].drop(labels=[s])
        if len(others) == 0 or others.isna().all():
            rows.append(dict(sample=s, nearest=None, correlation=None))
            continue
        nn = others.idxmax()
        rows.append(dict(sample=s, nearest=nn, correlation=round(float(others.max()), 4)))
    return pd.DataFrame(rows)


def flag_duplicates(corr: pd.DataFrame, thresh: float = 0.99) -> list:
    """Sample pairs correlated above `thresh` — likely relabels/duplicates."""
    dup = []
    cols = list(corr.columns)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            v = corr.iloc[i, j]
            if pd.notna(v) and v >= thresh:
                dup.append((cols[i], cols[j], round(float(v), 4)))
    return dup


def check_expected_pairs(corr: pd.DataFrame, expected_pairs: list, thresh: float = 0.99) -> list:
    """For known-identical pairs (e.g. a replicate), verify they are the mutual
    nearest neighbour and correlated above `thresh`; return any that fail."""
    nn = nearest_neighbours(corr).set_index("sample")
    problems = []
    for a, b in expected_pairs:
        if a not in corr.columns or b not in corr.columns:
            problems.append(dict(pair=(a, b), issue="sample not in matrix"))
            continue
        v = corr.at[a, b]
        mutual = (nn.at[a, "nearest"] == b) and (nn.at[b, "nearest"] == a)
        if pd.isna(v) or v < thresh or not mutual:
            problems.append(dict(pair=(a, b), correlation=None if pd.isna(v) else round(float(v), 4),
                                 mutual_nearest=bool(mutual), issue="pair not confidently identical"))
    return problems


def sex_check(sample_paths: dict) -> pd.DataFrame:
    """Per-sample methylation-based sex (chrY/autosome coverage ratio), reading
    each per-cytosine file. sample_paths: {sample: path}."""
    rows = []
    for s, p in sample_paths.items():
        try:
            meth = eio.read_methylation(p)
            cpg = meth[meth.context == "CpG"]
            rows.append(dict(sample=s, inferred_sex=eqc._infer_sex(cpg)))
        except Exception as e:
            rows.append(dict(sample=s, inferred_sex=f"error: {type(e).__name__}"))
    return pd.DataFrame(rows)


def summarise(corr: pd.DataFrame, expected_pairs: list | None = None) -> dict:
    dups = flag_duplicates(corr)
    problems = check_expected_pairs(corr, expected_pairs) if expected_pairs else []
    ok = (not problems) and (not dups or bool(expected_pairs))
    parts = []
    if dups:
        parts.append(f"{len(dups)} high-correlation pair(s) (>=0.99)")
    if problems:
        parts.append(f"{len(problems)} expected-identical pair(s) failed")
    return dict(n_samples=int(corr.shape[0]), duplicates=dups, expected_pair_problems=problems,
                ok=bool(ok and not problems),
                headline=("Identity checks passed." if not parts
                          else "Review: " + "; ".join(parts) + "."))
