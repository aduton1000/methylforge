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


def cmd_build_reference(a):
    from . import episignature as eepi
    mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
    lab = pd.read_csv(a.labels, sep="\t")
    labels = dict(zip(lab.iloc[:, 0], lab.iloc[:, 1]))
    b = eepi.build_reference(mat, labels, control_class=a.control_class,
                             per_class_n=a.features_per_class, min_delta=a.min_delta,
                             alpha=a.alpha, validate=not a.no_validate)
    out = eepi.save_reference(b, a.outdir)
    print(f"[methylforge] reference built: {len(b['features'])} features · classes "
          f"{sorted(set(labels.values()))} · saved to {a.outdir}")
    print(f"[methylforge] drop-in: --episign-ref {out['reference_matrix']} "
          f"--episign-labels {out['labels']}")
    if b["metrics"]:
        m = b["metrics"]
        print(f"[methylforge] leave-one-out CV accuracy {m['accuracy']} (n={m['n']}):")
        for c, pc in m["per_class"].items():
            print(f"    {c:18s} sensitivity {pc['sensitivity']}  specificity {pc['specificity']}  (n={pc['n']})")


def cmd_array(a):
    from . import array as earr
    r = earr.ingest(a.beta, manifest_path=a.manifest, out_matrix=a.out_matrix,
                    out_methyl_dir=a.out_methyl_dir, cov=a.cov)
    print(f"[methylforge] EPIC array ingested: {r['n_features']} features x {r['n_samples']} samples")
    if r.get("matrix"):
        print(f"[methylforge] coordinate-keyed matrix -> {r['matrix']}")
    if r.get("methyl_files"):
        print(f"[methylforge] per-sample methylation -> {a.out_methyl_dir} ({len(r['methyl_files'])} files)")


def cmd_cohort(a):
    from . import cohort as ecoh
    if a.samplesheet:
        paths, groups = ecoh.read_samplesheet(a.samplesheet)
    elif a.inputs:
        paths, groups = ecoh.paths_from_globs(a.inputs), {}
    else:
        raise SystemExit("provide --samplesheet or --inputs")
    sites = ecoh._load_sites(a.sites) if a.sites else None
    mat = ecoh.build_matrix(paths, min_cov=a.min_cov, min_samples_frac=a.min_samples_frac,
                            sites=sites, region_bed=a.region_bed, context=a.context)
    mat.to_csv(a.out, sep="\t")
    print(f"[methylforge] matrix {mat.shape[0]} features x {mat.shape[1]} samples -> {a.out}")
    if groups and a.groups_out:
        ecoh.write_groups(groups, a.groups_out)
        print(f"[methylforge] groups -> {a.groups_out}")


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

    br = sub.add_parser("build-reference", help="build + validate an episignature reference from a labelled cohort matrix")
    br.add_argument("--matrix", required=True, help="feature x sample beta matrix")
    br.add_argument("--labels", required=True, help="sample,class TSV (a 'Control' class required)")
    br.add_argument("-o", "--outdir", required=True)
    br.add_argument("--control-class", default="Control")
    br.add_argument("--features-per-class", type=int, default=40)
    br.add_argument("--min-delta", type=float, default=0.15)
    br.add_argument("--alpha", type=float, default=0.05)
    br.add_argument("--no-validate", action="store_true", help="skip leave-one-out CV")
    br.set_defaults(func=cmd_build_reference)

    ar = sub.add_parser("array", help="ingest an Illumina EPIC/450K beta matrix -> coord-keyed matrix / per-sample methyl")
    ar.add_argument("--beta", required=True, help="probe x sample beta matrix (from minfi/sesame)")
    ar.add_argument("--manifest", help="probe->chrom,pos TSV (omit if beta already coord-keyed)")
    ar.add_argument("--out-matrix", help="write coordinate-keyed feature x sample matrix")
    ar.add_argument("--out-methyl-dir", help="write per-sample methylation TSVs (for run/imprinting)")
    ar.add_argument("--cov", type=int, default=30, help="nominal coverage assigned to array betas")
    ar.set_defaults(func=cmd_array)

    co = sub.add_parser("cohort", help="build a feature x sample methylation matrix from many samples")
    co.add_argument("--samplesheet", help="CSV: sample,methylation[,group]")
    co.add_argument("--inputs", nargs="*", help="methylation files/globs (sample = basename)")
    co.add_argument("--out", required=True)
    co.add_argument("--groups-out", help="write sample,group TSV (needs group column)")
    co.add_argument("--sites", help="restrict to a feature list (chrom:pos, or chrom<TAB>pos)")
    co.add_argument("--region-bed", help="restrict to regions (targeted panel)")
    co.add_argument("--min-cov", type=int, default=5)
    co.add_argument("--min-samples-frac", type=float, default=0.75)
    co.add_argument("--context", default="CpG")
    co.set_defaults(func=cmd_cohort)

    t = sub.add_parser("testdata"); t.add_argument("-o", "--outdir", default="methylforge_testdata")
    t.set_defaults(func=cmd_testdata)

    v = sub.add_parser("version"); v.set_defaults(func=cmd_version)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
