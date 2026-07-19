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
