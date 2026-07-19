"""Stage 7 — episignature classification.

Reference-based classification of a genome-wide methylation profile against a
database of disorder-specific episignatures (following the framework of
Aref-Eshghi/Sadikovic). Each disorder's signature is a set of feature CpGs with
class-specific methylation. A query sample is scored against every class
centroid; the call is the best-matching class, with a calibrated confidence and
an explicit 'inconclusive' zone so uncertain samples are never over-called.

The default classifier is a nearest-centroid / correlation model (no sklearn
dependency, so it runs anywhere). A scikit-learn SVM back-end is used
automatically if scikit-learn is installed (production HPC env).

Reference format (TSV): first column 'feature' (CpG id), one column per
reference sample; plus a companion label map {sample: class}. Or a pre-computed
centroid table (feature x class).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

CONF_CALL = 0.80      # top class probability to make a positive call
CONF_MARGIN = 0.15    # min margin over 2nd class


def build_centroids(ref_matrix: pd.DataFrame, labels: dict) -> pd.DataFrame:
    """ref_matrix: feature x sample beta; labels: {sample: class} -> feature x class centroid."""
    classes = sorted(set(labels.values()))
    cols = {}
    for c in classes:
        samples = [s for s, lab in labels.items() if lab == c and s in ref_matrix.columns]
        cols[c] = ref_matrix[samples].mean(axis=1)
    return pd.DataFrame(cols)


def _try_sklearn(ref_matrix, labels, query):
    try:
        from sklearn.svm import SVC
    except Exception:
        return None
    samples = [s for s in ref_matrix.columns if s in labels]
    X = ref_matrix[samples].T.to_numpy(float)
    y = [labels[s] for s in samples]
    clf = SVC(kernel="linear", probability=True).fit(X, y)
    q = query.reindex(ref_matrix.index).to_numpy(float).reshape(1, -1)
    proba = clf.predict_proba(q)[0]
    return dict(zip(clf.classes_, proba))


def classify(query: pd.Series, centroids: pd.DataFrame,
             ref_matrix: pd.DataFrame | None = None, labels: dict | None = None) -> dict:
    """Classify a query methylation profile. query indexed by feature (beta 0..1)."""
    feats = centroids.index
    q = query.reindex(feats).to_numpy(float)
    mask = ~np.isnan(q)
    q_use = q[mask]

    probs = None
    method = "nearest-centroid (correlation)"
    if ref_matrix is not None and labels is not None:
        probs = _try_sklearn(ref_matrix, labels, query)
        if probs is not None:
            method = "scikit-learn linear SVM"

    if probs is None:
        # correlation of query to each class centroid over shared features
        corr = {}
        for c in centroids.columns:
            v = centroids[c].to_numpy(float)[mask]
            if v.std() > 0 and q_use.std() > 0:
                corr[c] = float(np.corrcoef(q_use, v)[0, 1])
            else:
                corr[c] = -1.0
        # softmax over correlations (temperature scaled) -> pseudo-probabilities
        z = np.array([corr[c] for c in centroids.columns]) / 0.10
        z = z - z.max()
        p = np.exp(z) / np.exp(z).sum()
        probs = dict(zip(centroids.columns, p))
        extra = {"correlation": {k: round(v, 3) for k, v in corr.items()}}
    else:
        extra = {}

    ranked = sorted(probs.items(), key=lambda kv: kv[1], reverse=True)
    top, top_p = ranked[0]
    second_p = ranked[1][1] if len(ranked) > 1 else 0.0
    margin = top_p - second_p

    called = (top_p >= CONF_CALL and margin >= CONF_MARGIN)
    if not called:
        verdict = "inconclusive"
        interp = ("No episignature met the calling threshold; profile does not "
                  "confidently match a known signature (or is control-like).")
        call_class = "inconclusive"
    elif top.lower().startswith("control"):
        verdict = "control-like"
        interp = "Methylation profile matches the control class — no disorder episignature detected."
        call_class = top
    else:
        verdict = "positive"
        interp = (f"Methylation profile matches the {top} episignature "
                  f"(p={top_p:.2f}, margin {margin:.2f}) — supports the associated diagnosis "
                  f"and can provide functional (PS3-type) evidence for a related variant.")
        call_class = top

    return dict(call=call_class, verdict=verdict, top_class=top,
                probability=round(top_p, 3), margin=round(margin, 3),
                method=method, interpretation=interp,
                probabilities={k: round(v, 3) for k, v in probs.items()},
                n_features_used=int(mask.sum()), **extra)


# ─────────────────────────────────────────────────────────────────────────────
# Reference builder — feature selection + leave-one-out validation
# ─────────────────────────────────────────────────────────────────────────────
import os as _os
import json as _json
from scipy import stats as _stats
from . import dmr as _dmr


def select_features(matrix: pd.DataFrame, labels: dict, control_class="Control",
                    per_class_n=40, min_delta=0.15, alpha=0.05) -> tuple[dict, list]:
    """Per non-control class, select the most discriminative CpGs vs controls
    (|Δβ| >= min_delta AND M-value Welch-t p < alpha), top `per_class_n` by |Δβ|.
    Returns ({class: [features]}, combined_feature_list)."""
    ctrl = [s for s, l in labels.items() if l == control_class and s in matrix.columns]
    if not ctrl:
        raise ValueError(f"no samples of control class '{control_class}'")
    K = matrix[ctrl].to_numpy(float)
    sig = {}
    for cls in sorted(set(labels.values())):
        if cls == control_class:
            continue
        cs = [s for s, l in labels.items() if l == cls and s in matrix.columns]
        C = matrix[cs].to_numpy(float)
        delta = np.nanmean(C, axis=1) - np.nanmean(K, axis=1)
        p = np.ones(len(matrix))
        for i in range(len(matrix)):
            a, b = C[i][~np.isnan(C[i])], K[i][~np.isnan(K[i])]
            if len(a) >= 2 and len(b) >= 2:
                try:
                    p[i] = _stats.ttest_ind(_dmr._mvalue(a), _dmr._mvalue(b), equal_var=False).pvalue
                except Exception:
                    p[i] = 1.0
        mask = (np.abs(delta) >= min_delta) & (p < alpha)
        order = np.argsort(-np.abs(delta * mask))
        feats = [matrix.index[i] for i in order if mask[i]][:per_class_n]
        sig[cls] = feats
    combined = sorted({f for fs in sig.values() for f in fs})
    return sig, combined


def _metrics(y_true, y_pred, classes):
    labs = sorted(set(classes))
    cols = labs + ["inconclusive"]
    conf = {t: {p: 0 for p in cols} for t in labs}
    for t, p in zip(y_true, y_pred):
        conf[t][p if p in cols else "inconclusive"] += 1
    per = {}
    for c in labs:
        tp = conf[c][c]
        fn = sum(v for p, v in conf[c].items() if p != c)
        fp = sum(conf[t][c] for t in labs if t != c)
        tn = sum(conf[t][p] for t in labs for p in cols if t != c and p != c)
        per[c] = dict(n=tp + fn,
                      sensitivity=round(tp / (tp + fn), 3) if (tp + fn) else None,
                      specificity=round(tn / (tn + fp), 3) if (tn + fp) else None,
                      tp=tp, fn=fn, fp=fp)
    acc = sum(1 for t, p in zip(y_true, y_pred) if t == p) / max(len(y_true), 1)
    return dict(accuracy=round(acc, 3), n=len(y_true), classes=labs,
                per_class=per, confusion=conf)


def cross_validate(matrix: pd.DataFrame, labels: dict, control_class="Control",
                   per_class_n=40, min_delta=0.15, alpha=0.05) -> dict:
    """Nested leave-one-out CV: features are re-selected on each training fold
    (no leakage), the held-out sample is classified, metrics are aggregated."""
    samples = [s for s in matrix.columns if s in labels]
    y_true, y_pred = [], []
    for held in samples:
        train = [s for s in samples if s != held]
        tl = {s: labels[s] for s in train}
        try:
            _, feats = select_features(matrix[train], tl, control_class, per_class_n, min_delta, alpha)
            if not feats:
                y_true.append(labels[held]); y_pred.append("inconclusive"); continue
            cent = build_centroids(matrix.loc[feats, train], tl)
            res = classify(matrix.loc[feats, held], cent, matrix.loc[feats, train], tl)
            y_pred.append(res["call"])
        except Exception:
            y_pred.append("inconclusive")
        y_true.append(labels[held])
    return _metrics(y_true, y_pred, [labels[s] for s in samples])


def build_reference(matrix: pd.DataFrame, labels: dict, control_class="Control",
                    per_class_n=40, min_delta=0.15, alpha=0.05, validate=True) -> dict:
    sig, feats = select_features(matrix, labels, control_class, per_class_n, min_delta, alpha)
    ref_matrix = matrix.loc[feats]
    centroids = build_centroids(ref_matrix, labels)
    metrics = cross_validate(matrix, labels, control_class, per_class_n, min_delta, alpha) if validate else None
    return dict(features=feats, signatures=sig, ref_matrix=ref_matrix,
                labels=labels, centroids=centroids, control_class=control_class,
                params=dict(per_class_n=per_class_n, min_delta=min_delta, alpha=alpha),
                metrics=metrics)


def save_reference(bundle: dict, outdir: str) -> dict:
    """Persist a reference bundle. reference_matrix.tsv + labels.tsv are drop-in
    for `methylforge run --episign-ref/--episign-labels`."""
    _os.makedirs(outdir, exist_ok=True)
    rm = _os.path.join(outdir, "reference_matrix.tsv")
    bundle["ref_matrix"].to_csv(rm, sep="\t")
    lab = _os.path.join(outdir, "labels.tsv")
    pd.DataFrame({"sample": list(bundle["labels"]), "class": list(bundle["labels"].values())}).to_csv(lab, sep="\t", index=False)
    bundle["centroids"].to_csv(_os.path.join(outdir, "centroids.tsv"), sep="\t")
    with open(_os.path.join(outdir, "signatures.json"), "w") as fh:
        _json.dump(bundle["signatures"], fh, indent=2)
    with open(_os.path.join(outdir, "manifest.json"), "w") as fh:
        _json.dump({"control_class": bundle["control_class"], "params": bundle["params"],
                    "n_features": len(bundle["features"]),
                    "classes": sorted(set(bundle["labels"].values())),
                    "metrics": bundle["metrics"]}, fh, indent=2)
    return dict(reference_matrix=rm, labels=lab,
                manifest=_os.path.join(outdir, "manifest.json"))
