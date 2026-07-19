"""Synthetic test-data generator for MethylForge.

Deterministic (fixed RNG seed) synthetic methylation data that exercises every
analysis stage and produces the *expected* biological result, so a fresh install
can be validated end-to-end without real sequencing data:

  Genome-wide per-CpG methylation (Bismark .cov) for three samples:
    * normal_control  — all imprinted loci normal, control episignature
    * pws_case        — SNRPN hypermethylated  -> expected: Prader-Willi call
    * sotos_case      — imprinted loci normal, Sotos_NSD1 episignature positive
  Each sample also contains the episignature feature CpGs and CHH sites (for the
  conversion-rate estimate) and chrY CpGs (for sex inference).

  Episignature reference matrix + labels (4 classes) and matching queries.
  A cohort DMR matrix (cases vs controls) with one implanted hypermethylated DMR.
  A small genes BED for DMR annotation.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from . import io as eio

RNG = np.random.default_rng(20260719)
EPI_CLASSES = ["Control", "Sotos_NSD1", "Kabuki_KMT2D", "CHARGE_CHD7"]
N_EPI_FEATURES = 60


def _cov_rows(records):
    """records: list of (chrom,start,beta,cov,context) -> bismark .cov dataframe (CpG lines only
    written to .cov; CHH kept separately handled by writing a native tsv instead)."""
    rows = []
    for chrom, start, beta, cov, ctx in records:
        meth = int(round(beta * cov))
        rows.append((chrom, start, start + 1, round(beta * 100, 3), meth, cov - meth, ctx))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "pct", "meth", "unmeth", "context"])


def _epi_feature_ids():
    return [f"chr{((i % 22) + 1)}:{100000 + i * 500}" for i in range(N_EPI_FEATURES)]


def _episig_reference():
    feats = _epi_feature_ids()
    # class-specific mean methylation over features (control ~0.5 baseline)
    base = RNG.uniform(0.4, 0.6, size=len(feats))
    centroids = {}
    for ci, cls in enumerate(EPI_CLASSES):
        v = base.copy()
        if cls != "Control":
            # each disorder perturbs a distinct block of features
            block = slice(ci * 12, ci * 12 + 12)
            v[block] = np.clip(base[block] + (0.35 if ci % 2 else -0.35), 0.02, 0.98)
        centroids[cls] = v
    # reference samples: 6 per class with noise
    cols = {}
    labels = {}
    for cls in EPI_CLASSES:
        for k in range(6):
            s = f"{cls}_ref{k+1}"
            cols[s] = np.clip(centroids[cls] + RNG.normal(0, 0.03, len(feats)), 0, 1)
            labels[s] = cls
    ref = pd.DataFrame(cols, index=feats)
    return ref, labels, centroids, feats


def _sample_records(imprint_overrides: dict, epi_class: str, ref_centroids, feats,
                    sex="female"):
    """Build a full per-cytosine record list for one sample."""
    panel = eio.read_dmr_panel()
    rec = []
    # imprinted DMRs: ~8 CpGs each, normal ~0.5 unless overridden
    for _, r in panel.iterrows():
        beta = imprint_overrides.get(r.locus, 0.5)
        for j in range(8):
            pos = int(r.start) + 40 * j
            rec.append((r.chrom, pos, float(np.clip(beta + RNG.normal(0, 0.03), 0, 1)), 30, "CpG"))
    # episignature feature CpGs (match chosen class centroid)
    cvec = ref_centroids[epi_class]
    for f, b in zip(feats, cvec):
        chrom, pos = f.split(":")
        rec.append((chrom, int(pos), float(np.clip(b + RNG.normal(0, 0.02), 0, 1)), 30, "CpG"))
    # background CpGs (bimodal genome-wide methylation)
    for _ in range(3000):
        chrom = f"chr{RNG.integers(1, 23)}"
        pos = int(RNG.integers(3_000_000, 3_500_000))
        b = RNG.choice([RNG.uniform(0, 0.15), RNG.uniform(0.8, 1.0)], p=[0.3, 0.7])
        rec.append((chrom, pos, float(b), int(RNG.integers(12, 45)), "CpG"))
    # CHH sites (near-zero methylation -> ~99% conversion)
    for _ in range(1500):
        chrom = f"chr{RNG.integers(1, 23)}"
        pos = int(RNG.integers(4_000_000, 4_200_000))
        rec.append((chrom, pos, float(np.clip(RNG.normal(0.006, 0.004), 0, 0.05)), 30, "CHH"))
    # chrY CpGs for sex inference
    ycov = 28 if sex == "male" else 1
    for _ in range(200):
        rec.append(("chrY", int(RNG.integers(2_700_000, 2_800_000)),
                    float(RNG.uniform(0, 1)), ycov, "CpG"))
    return rec


def _write_native(records, path):
    """Write full per-cytosine (incl. CHH) as MethylForge native TSV so QC can see CHH."""
    df = _cov_rows(records)
    df[["chrom", "start", "end", "context", "meth"]].assign(
        cov=df.meth + df.unmeth, beta=(df.meth / (df.meth + df.unmeth)).round(4)
    )[["chrom", "start", "end", "context", "meth", "cov", "beta"]].to_csv(path, sep="\t", index=False)


def generate(outdir: str) -> dict:
    os.makedirs(outdir, exist_ok=True)
    data = os.path.join(outdir, "data")
    os.makedirs(data, exist_ok=True)
    ref, labels, centroids, feats = _episig_reference()

    samples = {
        "normal_control": dict(overrides={}, epi="Control", sex="female"),
        "pws_case": dict(overrides={"PWS/AS-IC": 0.95}, epi="Control", sex="male"),
        "sotos_case": dict(overrides={}, epi="Sotos_NSD1", sex="female"),
    }
    paths = {}
    for name, spec in samples.items():
        rec = _sample_records(spec["overrides"], spec["epi"], centroids, feats, spec["sex"])
        p = os.path.join(data, f"{name}.methyl.tsv")
        _write_native(rec, p)
        paths[name] = p

    # episignature reference + labels + standalone queries
    ref.to_csv(os.path.join(data, "episig_reference.tsv"), sep="\t")
    pd.DataFrame({"sample": list(labels), "class": list(labels.values())}).to_csv(
        os.path.join(data, "episig_labels.tsv"), sep="\t", index=False)
    q_pos = pd.Series(np.clip(centroids["Sotos_NSD1"] + RNG.normal(0, 0.02, len(feats)), 0, 1),
                      index=feats, name="beta")
    q_pos.to_csv(os.path.join(data, "episig_query_sotos.tsv"), sep="\t")
    q_ctrl = pd.Series(np.clip(centroids["Control"] + RNG.normal(0, 0.02, len(feats)), 0, 1),
                       index=feats, name="beta")
    q_ctrl.to_csv(os.path.join(data, "episig_query_control.tsv"), sep="\t")

    # cohort DMR matrix: implant one hypermethylated DMR (chr3:5,000,000 region) in cases
    sites = [f"chr3:{5_000_000 + i * 60}" for i in range(12)] + \
            [f"chr3:{6_000_000 + i * 60}" for i in range(200)]
    cases = [f"case{i+1}" for i in range(6)]
    ctrls = [f"ctrl{i+1}" for i in range(6)]
    mat = pd.DataFrame(index=sites)
    for s in cases + ctrls:
        base = np.clip(RNG.uniform(0.1, 0.9, len(sites)), 0, 1)
        mat[s] = base
    for s in cases:  # DMR: first 12 sites hypermethylated in cases
        mat.loc[sites[:12], s] = np.clip(0.85 + RNG.normal(0, 0.04, 12), 0, 1)
    for s in ctrls:
        mat.loc[sites[:12], s] = np.clip(0.30 + RNG.normal(0, 0.04, 12), 0, 1)
    mat.round(4).to_csv(os.path.join(data, "dmr_matrix.tsv"), sep="\t")
    pd.DataFrame({"sample": cases + ctrls,
                  "group": ["case"] * 6 + ["control"] * 6}).to_csv(
        os.path.join(data, "dmr_groups.tsv"), sep="\t", index=False)

    # genes BED for DMR annotation
    with open(os.path.join(data, "genes.bed"), "w") as fh:
        fh.write("chr3\t4999000\t5002000\tTESTGENE1\n")
        fh.write("chr3\t5999000\t6002000\tTESTGENE2\n")

    print(f"[methylforge] test data written to {data}")
    for n, p in paths.items():
        print(f"  sample {n:16s} -> {p}")
    return dict(samples=paths, data=data,
                episig_reference=os.path.join(data, "episig_reference.tsv"),
                episig_labels=os.path.join(data, "episig_labels.tsv"),
                dmr_matrix=os.path.join(data, "dmr_matrix.tsv"),
                dmr_groups=os.path.join(data, "dmr_groups.tsv"),
                genes_bed=os.path.join(data, "genes.bed"))


if __name__ == "__main__":
    import sys
    generate(sys.argv[1] if len(sys.argv) > 1 else "methylforge_testdata")
