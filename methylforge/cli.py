"""MethylForge command-line interface.

Subcommands:
  run         end-to-end analysis of one sample -> clinical + research HTML,
              results.json, and a FamilyForge-compatible findings export
  qc          methylation QC only
  imprinting  imprinted-DMR analysis only
  dmr         differential methylation (cohort matrix)
  episign     episignature classification (reference + query)
  testdata    generate synthetic test data
  version
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import pandas as pd
from . import __version__, io as eio, qc as eqc, imprinting as eimp, \
    dmr as edmr, episignature as eepi, annotate as eann, pipeline as epipe
from . import export_familyforge as effx


def _p(obj):
    print(json.dumps(obj, indent=2, default=str))


def cmd_run(a):
    res = epipe.run(a.methylation, a.sample, a.outdir, family=a.family, assay=a.assay, genome=a.genome,
                    panel_path=a.panel, episign_ref=a.episign_ref,
                    episign_labels=a.episign_labels, dmr_matrix=a.dmr_matrix,
                    dmr_groups=a.dmr_groups, genes_bed=a.genes_bed)
    ff = os.path.join(a.outdir, f"{a.sample}.familyforge.json")
    effx.export(res, ff, family_id=a.family or a.sample)
    print(f"[methylforge] reports : {', '.join(res['outputs'].values())}")
    print(f"[methylforge] results : {os.path.join(a.outdir, a.sample + '.results.json')}")
    print(f"[methylforge] familyforge findings : {ff}")


def cmd_qc(a):
    _p(eqc.methylation_qc(eio.read_methylation(a.methylation)))


def cmd_imprinting(a):
    meth = eio.read_methylation(a.methylation)
    res = eimp.analyse_imprinting(meth, eio.read_dmr_panel(a.panel))
    print(res.to_string(index=False))
    _p(eimp.summarise_imprinting(res))


def cmd_dmr(a):
    mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
    gdf = pd.read_csv(a.groups, sep="\t")
    groups = dict(zip(gdf.iloc[:, 0], gdf.iloc[:, 1]))
    dmps = edmr.call_dmps(mat, groups)
    dmrs = edmr.merge_dmrs(dmps)
    if a.genes_bed:
        dmrs = eann.annotate_regions(dmrs, eann.load_genes_bed(a.genes_bed))
    _p(edmr.summarise_dmr(dmps, dmrs))
    if len(dmrs):
        print(dmrs.to_string(index=False))


def cmd_episign(a):
    ref = pd.read_csv(a.reference, sep="\t", index_col=0)
    lab = pd.read_csv(a.labels, sep="\t")
    labels = dict(zip(lab.iloc[:, 0], lab.iloc[:, 1]))
    centroids = eepi.build_centroids(ref, labels)
    q = pd.read_csv(a.query, sep="\t", index_col=0).iloc[:, 0]
    _p(eepi.classify(q, centroids, ref, labels))


def cmd_testdata(a):
    from . import testdata
    testdata.generate(a.outdir)


def cmd_version(a):
    print(f"MethylForge v{__version__}")


def build_parser():
    p = argparse.ArgumentParser(prog="methylforge",
                                description="MethylForge — epigenome (DNA-methylation) analysis pipeline")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="end-to-end single-sample analysis")
    r.add_argument("methylation", help="per-cytosine methylation file (.cov/.bedmethyl/.bedgraph/.tsv)")
    r.add_argument("-s", "--sample", required=True)
    r.add_argument("-o", "--outdir", default="methylforge_out")
    r.add_argument("--family", help="family id for FamilyForge export (default: sample)")
    r.add_argument("--assay", default="methylation")
    r.add_argument("--genome", default="GRCh38")
    r.add_argument("--panel", help="imprinted-DMR panel TSV (default: bundled)")
    r.add_argument("--episign-ref"); r.add_argument("--episign-labels")
    r.add_argument("--dmr-matrix"); r.add_argument("--dmr-groups")
    r.add_argument("--genes-bed")
    r.set_defaults(func=cmd_run)

    q = sub.add_parser("qc"); q.add_argument("methylation"); q.set_defaults(func=cmd_qc)

    im = sub.add_parser("imprinting"); im.add_argument("methylation")
    im.add_argument("--panel"); im.set_defaults(func=cmd_imprinting)

    d = sub.add_parser("dmr"); d.add_argument("--matrix", required=True)
    d.add_argument("--groups", required=True); d.add_argument("--genes-bed")
    d.set_defaults(func=cmd_dmr)

    e = sub.add_parser("episign"); e.add_argument("--reference", required=True)
    e.add_argument("--labels", required=True); e.add_argument("--query", required=True)
    e.set_defaults(func=cmd_episign)

    t = sub.add_parser("testdata"); t.add_argument("-o", "--outdir", default="methylforge_testdata")
    t.set_defaults(func=cmd_testdata)

    v = sub.add_parser("version"); v.set_defaults(func=cmd_version)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
