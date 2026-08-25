"""MethylForge self-test — generates synthetic data and asserts every analysis
stage produces the expected biological result. Run: `python tests/selftest.py`
or `pytest tests/selftest.py`."""
import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from methylforge import testdata, io as eio, qc as eqc, imprinting as eimp, \
    dmr as edmr, episignature as eepi, pipeline as epipe, export_familyforge as effx, \
    cohort as ecoh, array as earr, deconvolve as edc, clocks as eclk, asm as easm, \
    variability as evar, identity as eid


def test_all():
    tmp = tempfile.mkdtemp(prefix="methylforge_selftest_")
    d = testdata.generate(tmp)

    # QC
    meth = eio.read_methylation(d["samples"]["pws_case"])
    m = eqc.methylation_qc(meth)
    assert m["qc_pass"], "QC should pass on clean synthetic data"
    assert m["conversion_efficiency"] >= 0.98, "conversion ~100% expected"
    assert m["conversion_source"].startswith("lambda spike-in"), "spike-in should drive conversion QC"
    assert "puc19_methylation" in m["spikein"] and m["checks"].get("overconversion"), \
        "pUC19 over-conversion control + gate expected"

    # sex inference (chrX:autosome depth — robust to chrY mis-mapping)
    assert m["inferred_sex"].startswith("male"), f"pws_case should be male: {m['inferred_sex']}"
    sx = eqc.infer_sex(meth[meth.context == "CpG"])
    assert sx["sex"] == "male" and sx["confidence"] == "high", sx
    mN = eqc.methylation_qc(eio.read_methylation(d["samples"]["normal_control"]))
    assert mN["inferred_sex"].startswith("female"), f"normal_control should be female: {mN['inferred_sex']}"

    # imprinting — PWS case
    imp = eimp.analyse_imprinting(meth)
    snrpn = imp[imp.locus == "PWS/AS-IC"].iloc[0]
    assert snrpn["call"] == "hypermethylation", "SNRPN should be hypermethylated in PWS case"
    assert eimp.summarise_imprinting(imp)["disorders"] == ["Prader-Willi syndrome"]

    # imprinting — normal
    impN = eimp.analyse_imprinting(eio.read_methylation(d["samples"]["normal_control"]))
    assert (impN["call"] != "hypermethylation").all() and (impN["call"] != "hypomethylation").all()

    # episignature — Sotos positive, control negative
    ref = pd.read_csv(d["episig_reference"], sep="\t", index_col=0)
    lab = pd.read_csv(d["episig_labels"], sep="\t")
    labels = dict(zip(lab.iloc[:, 0], lab.iloc[:, 1]))
    cent = eepi.build_centroids(ref, labels)
    q_sotos = pd.read_csv(os.path.join(d["data"], "episig_query_sotos.tsv"), sep="\t", index_col=0).iloc[:, 0]
    q_ctrl = pd.read_csv(os.path.join(d["data"], "episig_query_control.tsv"), sep="\t", index_col=0).iloc[:, 0]
    assert eepi.classify(q_sotos, cent, ref, labels)["verdict"] == "positive"
    assert eepi.classify(q_ctrl, cent, ref, labels)["verdict"] == "control-like"

    # DMR
    mat = pd.read_csv(d["dmr_matrix"], sep="\t", index_col=0)
    gdf = pd.read_csv(d["dmr_groups"], sep="\t")
    groups = dict(zip(gdf.iloc[:, 0], gdf.iloc[:, 1]))
    dmrs = edmr.merge_dmrs(edmr.call_dmps(mat, groups))
    assert len(dmrs) >= 1 and (dmrs.direction == "hyper").any(), "expected a hypermethylated DMR"

    # cohort matrix builder -> DMR (should recover implanted DMR block + SNRPN)
    cpaths, cgroups = ecoh.read_samplesheet(d["cohort_samplesheet"])
    cmat = ecoh.build_matrix(cpaths)
    assert cmat.shape[1] == len(cpaths) and cmat.shape[0] > 50, "cohort matrix shape"
    cdmrs = edmr.merge_dmrs(edmr.call_dmps(cmat, cgroups))
    assert (cdmrs.direction == "hyper").sum() >= 2, "cohort DMR should recover implanted signals"

    # cohort QC dashboard — should flag the conversion-outlier control and render
    from methylforge import dashboard as edash
    qdf = edash.collect_qc(cpaths)
    assert (~qdf["qc_pass"]).sum() >= 1, "dashboard should flag the conversion outlier"
    dash_html = edash.render_dashboard(qdf, os.path.join(tmp, "cohort_qc.html"))
    assert os.path.getsize(dash_html) > 50000, "dashboard HTML should embed the overview figure"

    # episignature reference builder — leave-one-out CV on clean data should be strong
    bundle = eepi.build_reference(ref, labels, per_class_n=12, min_delta=0.15)
    assert bundle["metrics"]["accuracy"] >= 0.8, "episignature LOOCV accuracy should be high"
    assert len(bundle["features"]) > 0 and bundle["ref_matrix"].shape[1] == len(labels)

    # EPIC array ingestion -> coord matrix -> per-sample methyl -> imprinting
    coord = earr.map_probes_to_coords(earr.read_beta_matrix(d["epic_beta"]),
                                      earr.load_manifest(d["epic_manifest"]))
    assert coord.shape[1] == 3, "3 array samples expected"
    mp = os.path.join(tmp, "arr_pws.methyl.tsv")
    earr.to_sample_methyl(coord, "arr_pws", mp)
    imp_arr = eimp.analyse_imprinting(eio.read_methylation(mp))
    assert imp_arr[imp_arr.locus == "PWS/AS-IC"].iloc[0]["call"] == "hypermethylation", \
        "array PWS sample should flag SNRPN hypermethylation"

    # ── v0.2 native stages ──────────────────────────────────────────────────
    def _q(path):
        m = eio.read_methylation(path)
        return pd.Series({f"{c}:{s}": b for c, s, b in zip(m.chrom, m.start, m.beta)})

    # cell-type deconvolution — recover known mixture proportions (CP and RPC run)
    dref = edc.load_reference(d["cell_reference"])
    props = edc.deconvolve(_q(d["mixture"]), dref, method="CP")["proportions"]
    for ct, truth in d["true_proportions"].items():
        assert abs(props[ct] - truth) < 0.1, f"deconvolution off for {ct}: {props}"
    assert edc.deconvolve(_q(d["mixture"]), dref, method="RPC")["proportions"]  # RPC path runs

    # epigenetic clock — finite age, most features found
    ck = eclk.load_clock(d["clock"])
    age = eclk.predict_age(_q(d["samples"]["sotos_case"]), ck)
    assert age["dnam_age"] == age["dnam_age"] and age["frac_missing"] < 0.5, age
    assert len(eclk.cohort_acceleration([50, 60, 70], [45, 60, 75])) == 3

    # ASM — recover the implanted chr20 region, null region excluded
    hp1, hp2 = easm.load_pair(d["hp1"], d["hp2"])
    asites = easm.asm_from_haplotypes(hp1, hp2, min_cov=5)
    aregions = easm.merge_asm_regions(asites)
    assert len(aregions) >= 1 and (aregions.chrom == "chr20").any(), "expected a chr20 ASM region"
    assert (asites[asites.asm].pos < 5_500_000).all(), "ASM should be confined to the implanted region"
    assert (easm.candidate_asm_scan(hp1, min_cov=5, min_sites=4).chrom == "chr20").any()

    # differential variability (DVP) — runs and returns the schema
    dvps = evar.call_dvps(mat, groups)
    assert "dvp" in dvps.columns and "log2_var_ratio" in dvps.columns

    # covariate-adjusted & paired DMR — still recover the implanted DMR
    cov = pd.read_csv(d["dmr_covariates"], sep="\t").set_index("sample")
    pdf = pd.read_csv(d["dmr_pairs"], sep="\t"); pair = dict(zip(pdf.iloc[:, 0], pdf.iloc[:, 1]))
    dmrs_cov = edmr.merge_dmrs(edmr.call_dmps_glm(mat, groups, covariates=cov))
    assert (dmrs_cov.direction == "hyper").any(), "covariate DMR should recover the implanted DMR"
    assert len(edmr.merge_dmrs(edmr.call_dmps_glm(mat, groups, pair=pair))) >= 1, "paired DMR runs"

    # identity — detect an implanted duplicate column
    swapped = cmat.copy(); swapped["dupctrl"] = cmat["ctrl1"]
    dups = eid.flag_duplicates(eid.correlation_matrix(swapped))
    assert any("dupctrl" in (x, y) and "ctrl1" in (x, y) for x, y, _ in dups), "duplicate not detected"

    # config fix 1 — low-pass gates accept shallow coverage the default rejects
    shallow = pd.DataFrame({"chrom": ["chr1"] * 1200, "start": range(1200), "end": range(1, 1201),
                            "context": ["CpG"] * 1200, "meth": [3] * 1200, "cov": [5] * 1200,
                            "beta": [0.6] * 1200})
    assert not eqc.methylation_qc(shallow)["checks"]["coverage"], "5x should fail the standard gate"
    assert eqc.methylation_qc(shallow, gates=eqc.gates_for("low-pass"))["checks"]["coverage"], \
        "5x should pass the low-pass gate"
    # config fix 2 — autosomes-only drops chrX/Y
    tiny = pd.DataFrame({"case1": [0.9, 0.1], "ctrl1": [0.1, 0.1]}, index=["chr1:100", "chrX:200"])
    assert edmr.restrict_autosomes(tiny).index.tolist() == ["chr1:100"]

    # end-to-end + FamilyForge export
    res = epipe.run(d["samples"]["pws_case"], "selftest", os.path.join(tmp, "out"),
                    family="FAM", episign_ref=d["episig_reference"],
                    episign_labels=d["episig_labels"])
    for key, out in res["outputs"].items():
        if str(out).endswith(".html"):
            assert os.path.getsize(out) > 50000, "report HTML should be substantial (embedded figures)"
    # Excel findings workbook wired into the pipeline
    xlsx = res["outputs"].get("findings_xlsx")
    assert xlsx and os.path.isfile(xlsx), "run should emit a .findings.xlsx workbook"
    import openpyxl
    wb = openpyxl.load_workbook(xlsx)
    assert wb.sheetnames == ["Findings", "Imprinting loci", "QC"], "findings workbook sheets"
    ff = os.path.join(tmp, "out", "selftest.familyforge.json")
    effx.export(res, ff, family_id="FAM")
    import json
    payload = json.load(open(ff))
    assert any(r["model"] == "epi_imprinting" and r["tier"] == "A" for r in payload["rows"])
    print("MethylForge self-test PASSED — all stages produce expected results.")


if __name__ == "__main__":
    test_all()
