"""Report figures for MethylForge (matplotlib -> base64 PNG for HTML embedding)."""
from __future__ import annotations
import io as _io
import base64
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito
BLUE = "#0072B2"; ORANGE = "#E69F00"; GREEN = "#009E73"; RED = "#D55E00"
SKY = "#56B4E9"; PUR = "#CC79A7"; GREY = "#999999"; YEL = "#F0E442"; DARK = "#222222"
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})


def _b64(fig) -> str:
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def imprinting_plot(res) -> str:
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    y = np.arange(len(res))[::-1]
    colors = {"normal": GREEN, "hypermethylation": RED, "hypomethylation": BLUE,
              "no_call": GREY}
    betas = [b if b is not None else np.nan for b in res.observed_beta]
    ax.axvspan(0.35, 0.65, color=GREEN, alpha=0.10, label="normal range")
    ax.axvline(0.5, color=GREY, ls=":", lw=1)
    ax.barh(y, betas, color=[colors.get(c, GREY) for c in res.call], height=0.6,
            edgecolor=DARK, linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}\n({g})" for l, g in zip(res.locus, res.gene)], fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Methylation (beta)")
    ax.set_title("Imprinted-DMR methylation", fontweight="bold")
    for yi, (b, c) in enumerate(zip(betas[::-1], res.call[::-1])):
        if b == b:
            ax.text(min(b + 0.02, 0.96), y[::-1][yi], f"{b:.2f}", va="center", fontsize=7.5)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    return _b64(fig)


def dmr_volcano(dmps) -> str:
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    q = dmps.qvalue.clip(lower=1e-12)
    y = -np.log10(q)
    ax.scatter(dmps.delta, y, s=6, c=np.where(dmps.dmp, RED, GREY),
               alpha=0.6, edgecolors="none")
    ax.axhline(-np.log10(0.05), color=BLUE, ls="--", lw=1)
    ax.axvline(0.10, color=GREY, ls=":", lw=1); ax.axvline(-0.10, color=GREY, ls=":", lw=1)
    ax.set_xlabel("Methylation difference (case − control)")
    ax.set_ylabel("−log10(FDR q)")
    ax.set_title("Differential methylation (volcano)", fontweight="bold")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    return _b64(fig)


def episignature_plot(cls) -> str:
    probs = cls["probabilities"]
    items = sorted(probs.items(), key=lambda kv: kv[1], reverse=True)
    labels = [k for k, _ in items]
    vals = [v for _, v in items]
    cols = [RED if k == cls["top_class"] and cls["verdict"] == "positive"
            else (GREEN if k.lower().startswith("control") else GREY) for k in labels]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.bar(labels, vals, color=cols, edgecolor=DARK, linewidth=0.6)
    ax.axhline(0.80, color=BLUE, ls="--", lw=1, label="call threshold")
    ax.set_ylim(0, 1); ax.set_ylabel("Class probability")
    ax.set_title("Episignature classification", fontweight="bold")
    ax.tick_params(axis="x", rotation=30)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    return _b64(fig)


def qc_plot(meth, metrics) -> str:
    cpg = meth[(meth.context == "CpG") & (meth['cov'] > 0)]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    axes[0].hist(cpg['cov'].clip(upper=np.percentile(cpg['cov'], 99) if len(cpg) else 1),
                 bins=40, color=BLUE, edgecolor="white")
    axes[0].axvline(metrics["mean_cpg_coverage"], color=RED, ls="--",
                    label=f"mean {metrics['mean_cpg_coverage']:.0f}×")
    axes[0].set_xlabel("CpG coverage"); axes[0].set_ylabel("CpGs")
    axes[0].set_title("Coverage distribution", fontweight="bold"); axes[0].legend(frameon=False, fontsize=8)
    axes[1].hist(cpg.beta.dropna(), bins=40, color=GREEN, edgecolor="white")
    axes[1].set_xlabel("Methylation (beta)"); axes[1].set_ylabel("CpGs")
    axes[1].set_title("Global methylation", fontweight="bold")
    for ax in axes:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    return _b64(fig)


def cohort_qc_figure(df, gates) -> str:
    """Multi-panel batch overview (conversion, coverage, callable, global meth,
    spike-ins, pass/fail) across all samples in a run."""
    import pandas as _pd
    n = len(df)
    x = np.arange(n)
    labels = df["sample"].tolist()
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2))

    def bars(ax, vals, title, thr=None, thr_lab=None, good="high", ylabel=""):
        v = [0 if (x is None or (isinstance(x, float) and x != x)) else x for x in vals]
        if thr is None:
            cols = [BLUE] * n
        else:
            cols = [GREEN if ((val >= thr) if good == "high" else (val <= thr)) else RED for val in v]
        ax.bar(x, v, color=cols, edgecolor=DARK, linewidth=0.5)
        if thr is not None:
            ax.axhline(thr, color=GREY, ls="--", lw=1, label=thr_lab)
            ax.legend(fontsize=7, frameon=False)
        ax.set_title(title, fontweight="bold", fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel(ylabel, fontsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    bars(axes[0, 0], [(v * 100 if v is not None and v == v else None) for v in df["conversion"]],
         "Conversion efficiency (%)", gates["min_conversion"] * 100, f"≥{gates['min_conversion']*100:.0f}%", "high", "%")
    bars(axes[0, 1], df["mean_cov"], "Mean CpG coverage (×)", gates["min_mean_cov"], f"≥{gates['min_mean_cov']:.0f}×", "high", "×")
    bars(axes[0, 2], df["callable"], "Callable CpGs", gates["min_cpg_callable"], f"≥{gates['min_cpg_callable']}", "high")
    bars(axes[1, 0], [(v * 100 if v is not None and v == v else 0) for v in df["global_cpg_meth"]],
         "Global CpG methylation (%)", None, None, "high", "%")

    ax = axes[1, 1]
    lam = [(v * 100 if v is not None and v == v else 0) for v in df["lambda_conversion"]]
    puc = [(v * 100 if v is not None and v == v else 0) for v in df["puc19_overconversion"]]
    w = 0.4
    ax.bar(x - w / 2, lam, w, color=GREEN, label="λ conversion", edgecolor=DARK, linewidth=0.5)
    ax.bar(x + w / 2, puc, w, color=ORANGE, label="pUC19 over-conv.", edgecolor=DARK, linewidth=0.5)
    ax.set_title("Spike-in controls (%)", fontweight="bold", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.legend(fontsize=7, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = axes[1, 2]
    npass = int(df["qc_pass"].sum()); nrev = n - npass
    ax.bar(["PASS", "REVIEW"], [npass, nrev], color=[GREEN, RED], edgecolor=DARK)
    ax.set_title("QC outcome", fontweight="bold", fontsize=10)
    for i, val in enumerate([npass, nrev]):
        ax.text(i, val, str(val), ha="center", va="bottom", fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.tight_layout()
    return _b64(fig)
