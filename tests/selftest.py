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
    cohort as ecoh, array as earr


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
