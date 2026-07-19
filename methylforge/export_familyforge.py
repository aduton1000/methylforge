"""FamilyForge interoperability export.

Emits MethylForge findings in a form FamilyForge can consume for interpretation and
report generation:

  1. <sample>.familyforge.json  — a candidates.json-compatible object
     ({"source","sample","family","context_summary","rows","qc"}) whose `rows`
     carry the same base keys FamilyForge candidates use (gene, model, tier,
     tier_reason, confidence, segregation, zygosity, flags) plus epigenetic
     extras. FamilyForge's clinical_report can render these rows directly, or
     they can be merged into an existing candidates.json.

  2. <sample>.epi_findings.tsv  — a tolerant TSV (gene-keyed, like CallForge's
     cnv_calls.tsv) that FamilyForge's generic TSV loaders can ingest as an
     external-findings input.

Tiering mirrors FamilyForge's A/B/C convention:
  A = diagnostic-grade methylation defect at a known imprinting locus, or a
      positive disorder episignature (functional PS3-type evidence);
  B = imprinting/episignature finding needing confirmation, or a DMR over a
      disease gene;
  C = research DMRs and uncertain findings.
"""
from __future__ import annotations
import os
import json
import csv


def _imprinting_rows(res):
    rows = []
    imp = res["imprinting"]
    for r in imp.to_dict(orient="records"):
        if r["call"] in ("hypermethylation", "hypomethylation") and r["disorder"]:
            rows.append({
                "gene": r["gene"], "model": "epi_imprinting", "tier": "A",
                "tier_reason": (f"{r['call']} at {r['locus']} (beta {r['observed_beta']}, "
                                f"Δ {r['delta']}) — consistent with {r['disorder']}"),
                "confidence": "high", "segregation": "n/a (single-sample methylation)",
                "zygosity": "n/a", "inheritance_note": "imprinting/methylation defect",
                "flags": ["epigenetic", "imprinting"],
                "finding_type": "imprinting_dmr", "region": f"{r['chrom']}:{r['start']}-{r['end']}",
                "call": r["call"], "direction": ("hyper" if "hyper" in r["call"] else "hypo"),
                "disorder": r["disorder"], "metric": r["observed_beta"],
                "acmg_label": "Pathogenic methylation defect (locus-validated)",
                "source": "methylforge"})
    return rows


def _episign_rows(res):
    e = res.get("episign")
    if not e or e["verdict"] != "positive":
        return []
    gene = e["top_class"].split("_")[-1] if "_" in e["top_class"] else e["top_class"]
    return [{
        "gene": gene, "model": "epi_episignature", "tier": "A",
        "tier_reason": (f"positive {e['top_class']} episignature (p={e['probability']}, "
                        f"margin {e['margin']}) — supports diagnosis / PS3-type functional evidence"),
        "confidence": "high", "segregation": "n/a (methylation classifier)",
        "zygosity": "n/a", "inheritance_note": "episignature (genome-wide methylation)",
        "flags": ["epigenetic", "episignature", "PS3_functional"],
        "finding_type": "episignature", "region": "genome-wide",
        "call": "positive", "direction": "", "disorder": e["top_class"],
        "metric": e["probability"], "acmg_label": "Functional (PS3) support",
        "source": "methylforge"}]


def _dmr_rows(res, top=10):
    dmrs = res.get("dmrs")
    if dmrs is None or len(dmrs) == 0:
        return []
    rows = []
    for r in dmrs.head(top).to_dict(orient="records"):
        gene = r.get("nearest_gene", "") or "."
        rows.append({
            "gene": gene, "model": "epi_dmr", "tier": "C",
            "tier_reason": (f"{r['direction']}methylated DMR {r['chrom']}:{r['start']}-{r['end']} "
                            f"({r['n_cpg']} CpGs, mean Δ {r['mean_delta']}, q≤{r['min_qvalue']:.1e})"),
            "confidence": "moderate", "segregation": "n/a (cohort DMR)",
            "zygosity": "n/a", "inheritance_note": "differential methylation (research)",
            "flags": ["epigenetic", "dmr", "research"],
            "finding_type": "dmr", "region": f"{r['chrom']}:{r['start']}-{r['end']}",
            "call": r["direction"], "direction": r["direction"],
            "disorder": "", "metric": r["mean_delta"], "source": "methylforge"})
    return rows


def export(res: dict, out_json: str, family_id: str | None = None) -> str:
    rows = _imprinting_rows(res) + _episign_rows(res) + _dmr_rows(res)
    summary = {"n_candidates": len(rows),
               "tierA": sum(r["tier"] == "A" for r in rows),
               "tierB": sum(r["tier"] == "B" for r in rows),
               "tierC": sum(r["tier"] == "C" for r in rows)}
    payload = {
        "source": "methylforge", "methylforge_version": res.get("version", ""),
        "sample": res.get("sample"), "family": family_id or res.get("sample"),
        "assay": res.get("assay"), "genome": res.get("genome"),
        "provenance": res.get("provenance"),
        "context_summary": summary, "rows": rows,
        "imprinting_summary": res.get("imprint_summary"),
        "episignature": res.get("episign"),
        "qc": res.get("qc"),
    }
    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)

    # companion tolerant TSV
    tsv = os.path.splitext(out_json)[0].replace(".familyforge", "") + ".epi_findings.tsv"
    cols = ["sample", "family", "gene", "finding_type", "region", "call",
            "direction", "disorder", "tier", "metric", "tier_reason"]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({"sample": res.get("sample"), "family": family_id or res.get("sample"),
                        "gene": r["gene"], "finding_type": r["finding_type"],
                        "region": r["region"], "call": r["call"], "direction": r["direction"],
                        "disorder": r["disorder"], "tier": r["tier"], "metric": r["metric"],
                        "tier_reason": r["tier_reason"]})
    return out_json
