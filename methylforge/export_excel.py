"""Excel (.xlsx) findings export.

A spreadsheet listing of every MethylForge finding, for clinicians / geneticists
who work in Excel and for hand-off into LIMS. Two products, built from the same
tiered finding rows the FamilyForge export uses (so the Excel, the HTML report and
the FamilyForge hand-off never disagree):

  * per-sample   write_sample_findings_xlsx(res, path)
        sheets: Findings · Imprinting loci · QC
  * cohort       write_cohort_findings_xlsx(samples, path)
        sheets: All findings · Sample summary · QC

Findings are tier-coloured with the report's ACMG/AMP severity convention
(A pathogenic = red, B = amber, C/research = grey). Header uses GBCF navy.
"""
from __future__ import annotations
import os
import json
import datetime
import pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from . import export_familyforge as effx, __version__

NAVY = "16305E"; GOLD = "C9A24B"
TIER_FILL = {"A": "F4CCCC", "B": "FCE5CD", "C": "EDEDED"}   # pale red / amber / grey
TIER_FONT = {"A": "C0392B", "B": "9C6500", "C": "555555"}

_HEADER_FILL = PatternFill("solid", fgColor=NAVY)
_HEADER_FONT = Font(color="FFFFFF", bold=True, size=11)
_TITLE_FONT = Font(color=NAVY, bold=True, size=14)
_THIN = Side(style="thin", color="D0D0D0")
_BORDER = Border(left=_THIN, right=_THIN, top=_THIN, bottom=_THIN)

FINDING_COLS = ["sample", "family", "tier", "finding_type", "gene", "disorder",
                "call", "direction", "region", "metric", "confidence",
                "acmg_label", "tier_reason", "source"]


def _sample_finding_rows(res: dict) -> list[dict]:
    """All tiered findings for one sample (imprinting + episignature + DMR)."""
    rows = effx._imprinting_rows(res) + effx._episign_rows(res) + effx._dmr_rows(res)
    sample = res.get("sample"); family = res.get("family", sample)
    for r in rows:
        r.setdefault("acmg_label", "")
        r["sample"] = sample; r["family"] = family
    return rows


def _style_table(ws, df: pd.DataFrame, start_row: int, tier_col: str | None = "tier"):
    """Write a dataframe as a styled, auto-filtered, frozen table starting at start_row."""
    for j, col in enumerate(df.columns, start=1):
        c = ws.cell(row=start_row, column=j, value=str(col))
        c.fill = _HEADER_FILL; c.font = _HEADER_FONT; c.border = _BORDER
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for i, (_, rr) in enumerate(df.iterrows(), start=start_row + 1):
        tier = str(rr[tier_col]) if tier_col and tier_col in df.columns else None
        for j, col in enumerate(df.columns, start=1):
            val = rr[col]
            if isinstance(val, (list, tuple)):
                val = ", ".join(map(str, val))
            c = ws.cell(row=i, column=j, value=val)
            c.border = _BORDER
            c.alignment = Alignment(vertical="top", wrap_text=(col in ("tier_reason", "acmg_label", "interpretation")))
            if tier in TIER_FILL:
                c.fill = PatternFill("solid", fgColor=TIER_FILL[tier])
                if col == tier_col:
                    c.font = Font(bold=True, color=TIER_FONT[tier])
    n = len(df)
    ws.auto_filter.ref = f"A{start_row}:{get_column_letter(len(df.columns))}{start_row + n}"
    ws.freeze_panes = ws.cell(row=start_row + 1, column=1)
    # column widths
    for j, col in enumerate(df.columns, start=1):
        w = max([len(str(col))] + [len(str(v)) for v in df[col].astype(str).head(50)])
        ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 10), 60)


def _title(ws, text, sub=None):
    ws.cell(row=1, column=1, value=text).font = _TITLE_FONT
    if sub:
        ws.cell(row=2, column=1, value=sub).font = Font(color="666666", italic=True, size=9)
    return 4 if sub else 3


def _qc_kv(qc: dict) -> pd.DataFrame:
    spk = qc.get("spikein", {}) or {}
    kv = [
        ("QC pass", qc.get("qc_pass")),
        ("Conversion efficiency", qc.get("conversion_efficiency")),
        ("Conversion source", qc.get("conversion_source")),
        ("Mean CpG coverage", qc.get("mean_cpg_coverage")),
        ("Median CpG coverage", qc.get("median_cpg_coverage")),
        ("Callable CpGs", qc.get("n_cpg_callable")),
        ("Total CpGs", qc.get("n_cpg")),
        ("Global CpG methylation", qc.get("global_cpg_methylation")),
        ("Global CHH methylation", qc.get("global_chh_methylation")),
        ("Lambda conversion (spike-in)", spk.get("lambda_conversion")),
        ("pUC19 methylation (spike-in)", spk.get("puc19_methylation")),
        ("pUC19 over-conversion", spk.get("puc19_overconversion")),
        ("Inferred sex", qc.get("inferred_sex")),
    ]
    return pd.DataFrame(kv, columns=["metric", "value"])


def write_sample_findings_xlsx(res: dict, path: str) -> str:
    from openpyxl import Workbook
    wb = Workbook(); wb.remove(wb.active)
    sample = res.get("sample", "sample")

    # Findings
    ws = wb.create_sheet("Findings")
    r0 = _title(ws, f"MethylForge findings — {sample}",
                f"family {res.get('family', sample)} · {res.get('assay','methylation')} · "
                f"{res.get('genome','GRCh38')} · MethylForge v{__version__} · {datetime.date.today().isoformat()}")
    rows = _sample_finding_rows(res)
    if rows:
        df = pd.DataFrame(rows)[[c for c in FINDING_COLS if c in rows[0]]]
    else:
        df = pd.DataFrame([{**{c: "" for c in FINDING_COLS},
                            "sample": sample, "family": res.get("family", sample),
                            "tier": "", "finding_type": "none",
                            "tier_reason": "No reportable methylation findings at the loci examined."}])[FINDING_COLS]
    _style_table(ws, df, r0)

    # Imprinting loci (full panel)
    ws2 = wb.create_sheet("Imprinting loci")
    r0 = _title(ws2, f"Imprinted-DMR panel — {sample}")
    imp = res["imprinting"]
    idf = (imp if isinstance(imp, pd.DataFrame) else pd.DataFrame(imp))[
        ["locus", "gene", "chrom", "start", "end", "observed_beta", "normal_meth",
         "delta", "call", "disorder"]]
    _style_table(ws2, idf, r0, tier_col=None)

    # QC
    ws3 = wb.create_sheet("QC")
    r0 = _title(ws3, f"Methylation QC — {sample}")
    _style_table(ws3, _qc_kv(res["qc"]), r0, tier_col=None)

    wb.save(path)
    return path


def write_cohort_findings_xlsx(samples: list[dict], path: str,
                               title: str = "MethylForge cohort findings") -> str:
    """samples: list of FamilyForge-export dicts (each has sample, family, rows, qc,
    imprinting_summary, episignature). Builds All findings · Sample summary · QC."""
    from openpyxl import Workbook
    wb = Workbook(); wb.remove(wb.active)

    # All findings
    all_rows = []
    for s in samples:
        for r in s.get("rows", []):
            rr = {c: r.get(c, "") for c in FINDING_COLS}
            rr["sample"] = s.get("sample"); rr["family"] = s.get("family", s.get("sample"))
            all_rows.append(rr)
    ws = wb.create_sheet("All findings")
    r0 = _title(ws, title, f"{len(samples)} samples · {len(all_rows)} findings · "
                           f"MethylForge v{__version__} · {datetime.date.today().isoformat()}")
    fdf = pd.DataFrame(all_rows, columns=FINDING_COLS) if all_rows else \
        pd.DataFrame([{c: "" for c in FINDING_COLS}])
    if all_rows:
        fdf = fdf.sort_values(["tier", "family", "sample"]).reset_index(drop=True)
    _style_table(ws, fdf, r0)

    # Sample summary (one row per sample)
    summ = []
    for s in samples:
        cs = s.get("context_summary", {}) or {}
        isum = s.get("imprinting_summary", {}) or {}
        epi = s.get("episignature") or {}
        qc = s.get("qc", {}) or {}
        summ.append(dict(
            sample=s.get("sample"), family=s.get("family", s.get("sample")),
            headline=isum.get("headline", ""),
            imprinting_severity=isum.get("severity", ""),
            episignature=(f"{epi.get('verdict','')}"
                          + (f" ({epi.get('top_class')})" if epi.get("verdict") == "positive" else "")),
            tierA=cs.get("tierA", 0), tierB=cs.get("tierB", 0), tierC=cs.get("tierC", 0),
            n_findings=cs.get("n_candidates", 0),
            qc_pass=qc.get("qc_pass"), inferred_sex=(qc.get("inferred_sex", "") or "").split()[0]))
    ws2 = wb.create_sheet("Sample summary")
    r0 = _title(ws2, "Per-sample summary")
    _style_table(ws2, pd.DataFrame(summ), r0, tier_col=None)

    # QC
    qc_rows = []
    for s in samples:
        qc = s.get("qc", {}) or {}; spk = qc.get("spikein", {}) or {}
        qc_rows.append(dict(
            sample=s.get("sample"), qc_pass=qc.get("qc_pass"),
            conversion=qc.get("conversion_efficiency"), conversion_source=qc.get("conversion_source"),
            mean_cov=qc.get("mean_cpg_coverage"), callable=qc.get("n_cpg_callable"),
            global_cpg_meth=qc.get("global_cpg_methylation"),
            lambda_conversion=spk.get("lambda_conversion"),
            puc19_overconversion=spk.get("puc19_overconversion"),
            inferred_sex=(qc.get("inferred_sex", "") or "").split()[0]))
    ws3 = wb.create_sheet("QC")
    r0 = _title(ws3, "Per-sample QC")
    _style_table(ws3, pd.DataFrame(qc_rows), r0, tier_col=None)

    wb.save(path)
    return path


def cohort_from_familyforge(paths: list[str], out_xlsx: str,
                            title: str = "MethylForge cohort findings") -> str:
    samples = [json.load(open(p)) for p in paths]
    return write_cohort_findings_xlsx(samples, out_xlsx, title=title)
