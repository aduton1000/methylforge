"""MethylForge self-test — generates synthetic data and asserts every analysis
stage produces the expected biological result. Run: `python tests/selftest.py`
or `pytest tests/selftest.py`."""
import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from methylforge import testdata, io as eio, qc as eqc, imprinting as eimp, \
    dmr as edmr, episignature as eepi, pipeline as epipe, export_familyforge as effx


def test_all():
    tmp = tempfile.mkdtemp(prefix="methylforge_selftest_")
    d = testdata.generate(tmp)

    # QC
    meth = eio.read_methylation(d["samples"]["pws_case"])
    m = eqc.methylation_qc(meth)
    assert m["qc_pass"], "QC should pass on clean synthetic data"
    assert m["conversion_efficiency"] >= 0.98, "conversion ~100% expected (CHH near 0)"

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

    # end-to-end + FamilyForge export
    res = epipe.run(d["samples"]["pws_case"], "selftest", os.path.join(tmp, "out"),
                    family="FAM", episign_ref=d["episig_reference"],
                    episign_labels=d["episig_labels"])
    for html in res["outputs"].values():
        assert os.path.getsize(html) > 50000, "report HTML should be substantial (embedded figures)"
    ff = os.path.join(tmp, "out", "selftest.familyforge.json")
    effx.export(res, ff, family_id="FAM")
    import json
    payload = json.load(open(ff))
    assert any(r["model"] == "epi_imprinting" and r["tier"] == "A" for r in payload["rows"])
    print("MethylForge self-test PASSED — all stages produce expected results.")


if __name__ == "__main__":
    test_all()
