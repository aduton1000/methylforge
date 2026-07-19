"""Stage 9 — report engine.

Two report styles from one result set:
  * clinical  -> GBCF/Noguchi-branded clinical template (clinical_gbcf.html.j2)
                 with logos and the ACMG/AMP variant-severity colour convention
                 (Pathogenic=red, Likely pathogenic=orange->red, VUS=amber, Benign=green).
  * research  -> clean modern template (report.html.j2) with full DMR detail.
Also emits a machine-readable results.json.
"""
from __future__ import annotations
import os
import json
import base64
import datetime
from jinja2 import Environment, FileSystemLoader, select_autoescape
from . import plots, __version__

_TPL_DIR = os.path.join(os.path.dirname(__file__), "report_templates")
_ASSET_DIR = os.path.join(os.path.dirname(__file__), "resources", "assets")

_SEV_LABEL = {"p": "Pathogenic", "lp": "Likely pathogenic", "vus": "Uncertain (VUS)",
              "b": "Benign / normal", "none": "No call"}
_SEV_RANK = {"p": 4, "lp": 3, "vus": 2, "b": 1, "none": 0}


def _logo(name):
    p = os.path.join(_ASSET_DIR, name)
    if not os.path.isfile(p):
        return ""
    return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()


def _imprint_sev(call, disorder):
    if call in ("hypermethylation", "hypomethylation"):
        return ("p", _SEV_LABEL["p"]) if disorder else ("lp", _SEV_LABEL["lp"])
    if call == "normal":
        return ("b", _SEV_LABEL["b"])
    return ("none", _SEV_LABEL["none"])


def _episign_sev(verdict):
    return {"positive": "p", "control-like": "b", "inconclusive": "vus"}.get(verdict, "none")


def _assemble(results):
    """Attach severity classifications + build the clinical findings summary."""
    imp = results["imprinting"].copy()
    sevs, labels = [], []
    for _, r in imp.iterrows():
        s, l = _imprint_sev(r["call"], r["disorder"])
        sevs.append(s); labels.append(l)
    imp["_sev"] = sevs
    imp["_sev_label"] = labels

    episign = results.get("episign")
    findings = []
    for _, r in imp.iterrows():
        if r["call"] in ("hypermethylation", "hypomethylation"):
            findings.append(dict(type="Imprinting DMR",
                                 locus=f"{r['gene']} ({r['locus']})",
                                 obs=f"{r['call']} (β {r['observed_beta']})",
                                 disorder=r["disorder"], sev=r["_sev"], sev_label=r["_sev_label"]))
    if episign:
        es = _episign_sev(episign["verdict"])
        episign = {**episign, "_sev": f"sev-{es}"}
        if es in ("p", "vus"):
            findings.append(dict(type="Episignature", locus=episign["top_class"],
                                 obs=f"{episign['verdict']} (p={episign['probability']})",
                                 disorder=episign["top_class"] if es == "p" else "",
                                 sev=es, sev_label=_SEV_LABEL[es]))
    ranks = [_SEV_RANK[f["sev"]] for f in findings] or [0]
    top = max(ranks)
    headline_sev = "sev-" + next(k for k, v in _SEV_RANK.items() if v == top)
    return imp, episign, findings, headline_sev


def _headline(imprint_summary, episign):
    if imprint_summary["severity"] == "abnormal":
        return imprint_summary["headline"]
    if episign and episign["verdict"] == "positive":
        return episign["interpretation"]
    if episign and episign["verdict"] == "inconclusive":
        return "No abnormal imprinting methylation; episignature inconclusive — consider reflex testing."
    return "No abnormal methylation detected at the loci examined."


def render_report(results: dict, out_html: str, mode: str = "clinical",
                  date: str | None = None) -> str:
    env = Environment(loader=FileSystemLoader(_TPL_DIR),
                      autoescape=select_autoescape(["html", "xml"]))
    meth = results["meth"]
    qc = results["qc"]
    imp_summary = results["imprint_summary"]
    imp, episign, findings, headline_sev = _assemble(results)
    dmps = results.get("dmps")
    dmrs = results.get("dmrs")

    base = dict(
        sample=results.get("sample", "sample"),
        family=results.get("family", results.get("sample", "sample")),
        assay=results.get("assay", "methylation"),
        genome=results.get("genome", "GRCh38"),
        version=__version__,
        date=date or datetime.date.today().isoformat(),
        provenance=results.get("provenance", "methylforge run"),
        qc=qc, qc_plot=plots.qc_plot(meth, qc),
        imprinting=imp.to_dict(orient="records"),
        imprint_plot=plots.imprinting_plot(imp),
        episign=episign,
        episign_plot=plots.episignature_plot(episign) if episign else None,
    )

    if mode == "clinical":
        tpl = env.get_template("clinical_gbcf.html.j2")
        html = tpl.render(headline=_headline(imp_summary, episign), headline_sev=headline_sev,
                          findings=findings, logo_gbcf=_logo("gbcf_lab_logo.png"),
                          logo_noguchi=_logo("noguchi_logo.png"), **base)
    else:
        hc = {"sev-p": "abnormal", "sev-lp": "abnormal", "sev-vus": "inconclusive",
              "sev-b": "normal", "sev-none": "normal"}[headline_sev]
        tpl = env.get_template("report.html.j2")
        html = tpl.render(mode="research", headline=_headline(imp_summary, episign),
                          headline_class=hc, imprint_summary=imp_summary,
                          dmr_summary=results.get("dmr_summary"),
                          dmr_volcano=plots.dmr_volcano(dmps) if dmps is not None else None,
                          dmrs=dmrs.to_dict(orient="records") if dmrs is not None and len(dmrs) else [],
                          **base)
    with open(out_html, "w") as fh:
        fh.write(html)
    return out_html


def write_results_json(results: dict, out_json: str) -> str:
    def _clean(v):
        try:
            import pandas as pd
            if isinstance(v, pd.DataFrame):
                return v.to_dict(orient="records")
        except Exception:
            pass
        return v
    payload = {k: _clean(v) for k, v in results.items() if k != "meth"}
    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return out_json
