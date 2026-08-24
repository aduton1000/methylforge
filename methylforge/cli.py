"""MethylForge command-line interface.

Subcommands:
  run         end-to-end analysis of one sample -> clinical + research HTML,
              results.json, an Excel findings workbook, and a FamilyForge export
  findings-xlsx  aggregate per-sample FamilyForge findings into one Excel workbook
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
    clocks = None
    if getattr(a, "clock", None):
        clocks = {}
        for spec in a.clock:
            name, path = spec.split("=", 1) if "=" in spec else (os.path.basename(spec), spec)
            clocks[name] = path
    res = epipe.run(a.methylation, a.sample, a.outdir, family=a.family, assay=a.assay, genome=a.genome,
                    panel_path=a.panel, episign_ref=a.episign_ref,
                    episign_labels=a.episign_labels, dmr_matrix=a.dmr_matrix,
                    dmr_groups=a.dmr_groups, genes_bed=a.genes_bed,
                    gates_profile="low-pass" if getattr(a, "low_pass", False) else "standard",
                    cell_reference=getattr(a, "cell_reference", None), clocks=clocks,
                    chrono_age=getattr(a, "chrono_age", None),
                    dmr_covariates=getattr(a, "dmr_covariates", None),
                    dmr_pairs=getattr(a, "dmr_pairs", None),
                    dmr_autosomes=getattr(a, "dmr_autosomes", False),
                    asm_hp1=getattr(a, "asm_hp1", None), asm_hp2=getattr(a, "asm_hp2", None))
    ff = os.path.join(a.outdir, f"{a.sample}.familyforge.json")
    effx.export(res, ff, family_id=a.family or a.sample)
    print(f"[methylforge] reports : {', '.join(res['outputs'].values())}")
    print(f"[methylforge] results : {os.path.join(a.outdir, a.sample + '.results.json')}")
    print(f"[methylforge] familyforge findings : {ff}")


def cmd_qc(a):
    gates = eqc.gates_for("low-pass" if getattr(a, "low_pass", False) else "standard",
                          min_mean_cov=getattr(a, "min_mean_cov", None),
                          min_conversion=getattr(a, "min_conversion", None),
                          min_site_cov=getattr(a, "min_site_cov", None))
    _p(eqc.methylation_qc(eio.read_methylation(a.methylation), gates=gates))


def cmd_imprinting(a):
    meth = eio.read_methylation(a.methylation)
    res = eimp.analyse_imprinting(meth, eio.read_dmr_panel(a.panel))
    print(res.to_string(index=False))
    _p(eimp.summarise_imprinting(res))


def cmd_dmr(a):
    mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
    gdf = pd.read_csv(a.groups, sep="\t")
    groups = dict(zip(gdf.iloc[:, 0], gdf.iloc[:, 1]))
    if getattr(a, "autosomes", False):
        n0 = len(mat)
        mat = edmr.restrict_autosomes(mat)
        print(f"[methylforge] autosomes-only: {len(mat)}/{n0} CpGs retained (chrX/Y/M dropped)")
    covariates = pd.read_csv(a.covariates, sep=None, engine="python", index_col=0) \
        if getattr(a, "covariates", None) else None
    pair = None
    if getattr(a, "pair", None):
        pdf = pd.read_csv(a.pair, sep=None, engine="python")
        pair = dict(zip(pdf.iloc[:, 0], pdf.iloc[:, 1]))
    use_glm = covariates is not None or pair is not None or getattr(a, "model", "auto") == "glm"
    if use_glm:
        dmps = edmr.call_dmps_glm(mat, groups, covariates=covariates, pair=pair)
        model = "covariate/paired OLS on M-values"
    else:
        dmps = edmr.call_dmps(mat, groups)
        model = "two-group (Welch t / Mann-Whitney)"
    dmrs = edmr.merge_dmrs(dmps)
    if a.genes_bed:
        dmrs = eann.annotate_regions(dmrs, eann.load_genes_bed(a.genes_bed))
    print(f"[methylforge] DMR model: {model}")
    _p(edmr.summarise_dmr(dmps, dmrs))
    if len(dmrs):
        print(dmrs.to_string(index=False))


def _query_series(meth):
    """feature('chrom:pos') -> beta Series from a per-cytosine methylation frame."""
    q = pd.Series({f"{c}:{s}": b for c, s, b in zip(meth.chrom, meth.start, meth.beta)}, dtype=float)
    return q[~q.index.duplicated(keep="first")]


def cmd_deconvolve(a):
    from . import deconvolve as edc
    ref = edc.load_reference(a.reference)
    if a.matrix:
        mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
        props = edc.deconvolve_matrix(mat, ref, method=a.method)
        print(props.round(4).to_string())
        if a.out:
            props.round(4).to_csv(a.out, sep="\t")
            print(f"[methylforge] proportions -> {a.out}")
    else:
        res = edc.deconvolve(_query_series(eio.read_methylation(a.methylation)), ref, method=a.method)
        _p({**res, **edc.summarise(res)})


def cmd_clock(a):
    from . import clocks as eclk
    q = _query_series(eio.read_methylation(a.methylation))
    clocks = {}
    for spec in a.clock:
        name, path = spec.split("=", 1) if "=" in spec else (None, spec)
        ck = eclk.load_clock(path)
        if name:
            ck["meta"]["name"] = name
        clocks[ck["meta"]["name"]] = ck
    tbl = eclk.run_clocks(q, clocks)
    if a.chrono_age is not None:
        tbl["delta_age"] = [eclk.delta_age(x, a.chrono_age) for x in tbl["dnam_age"]]
    print(tbl.to_string(index=False))


def cmd_asm(a):
    from . import asm as easm
    if a.candidate:
        reg = easm.candidate_asm_scan(eio.read_methylation(a.candidate), min_cov=a.min_cov)
        print(f"[methylforge] candidate ASM (heuristic) regions: {len(reg)} "
              f"— confirm with phased/haplotype data")
        if len(reg):
            print(reg.to_string(index=False))
        return
    if a.allele_table:
        hp1, hp2, (n1, n2) = easm.read_allele_table(a.allele_table)
    elif a.hp1 and a.hp2:
        hp1, hp2 = easm.load_pair(a.hp1, a.hp2)
        n1, n2 = "hp1", "hp2"
    else:
        raise SystemExit("provide --hp1/--hp2, --allele-table, or --candidate")
    sites = easm.asm_from_haplotypes(hp1, hp2, min_cov=a.min_cov, min_delta=a.min_delta,
                                     fdr=a.fdr, hp1_name=n1, hp2_name=n2)
    regions = easm.merge_asm_regions(sites)
    _p(easm.summarise_asm(sites, regions))
    if len(regions):
        print(regions.to_string(index=False))
    if a.out:
        sites.to_csv(a.out, sep="\t", index=False)
        print(f"[methylforge] ASM sites -> {a.out}")


def cmd_dvp(a):
    from . import variability as evar
    mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
    gdf = pd.read_csv(a.groups, sep="\t")
    groups = dict(zip(gdf.iloc[:, 0], gdf.iloc[:, 1]))
    dvps = evar.call_dvps(mat, groups)
    _p(evar.summarise_dvp(dvps))
    top = dvps[dvps.dvp]
    if len(top):
        print(top.head(50).to_string(index=False))


def cmd_identity(a):
    from . import identity as eid, cohort as ecoh
    paths = None
    if a.matrix:
        mat = pd.read_csv(a.matrix, sep="\t", index_col=0)
    elif a.samplesheet:
        paths, _ = ecoh.read_samplesheet(a.samplesheet)
        mat = ecoh.build_matrix(paths)
    elif a.inputs:
        paths = ecoh.paths_from_globs(a.inputs)
        mat = ecoh.build_matrix(paths)
    else:
        raise SystemExit("provide --matrix, --samplesheet, or --inputs")
    corr = eid.correlation_matrix(mat)
    print("=== pairwise correlation ===")
    print(corr.round(4).to_string())
    print("=== nearest neighbours ===")
    print(eid.nearest_neighbours(corr).to_string(index=False))
    expected = [tuple(p.split(":")) for p in a.expected_pairs] if a.expected_pairs else None
    _p(eid.summarise(corr, expected))
    if a.sex_check and paths:
        print("=== methylation-based sex ===")
        print(eid.sex_check(paths).to_string(index=False))


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
                            sites=sites, region_bed=a.region_bed, context=a.context,
                            autosomes_only=getattr(a, "autosomes_only", False))
    mat.to_csv(a.out, sep="\t")
    print(f"[methylforge] matrix {mat.shape[0]} features x {mat.shape[1]} samples -> {a.out}")
    if groups and a.groups_out:
        ecoh.write_groups(groups, a.groups_out)
        print(f"[methylforge] groups -> {a.groups_out}")


def cmd_cohortqc(a):
    import glob
    from . import dashboard as edash, cohort as ecoh
    if a.samplesheet:
        paths, _ = ecoh.read_samplesheet(a.samplesheet)
        df = edash.collect_qc(paths)
    elif a.inputs:
        df = edash.collect_qc(ecoh.paths_from_globs(a.inputs))
    elif a.results:
        rp = [f for pat in a.results for f in sorted(glob.glob(pat))]
        df = edash.collect_from_results(rp)
    else:
        raise SystemExit("provide --samplesheet, --inputs, or --results")
    edash.render_dashboard(df, a.out, title=a.title)
    print(f"[methylforge] cohort QC dashboard ({len(df)} samples, "
          f"{int(df['qc_pass'].sum())} pass / {int((~df['qc_pass']).sum())} review) -> {a.out}")
    if a.summary_out:
        edash.write_summary(df, a.summary_out)
        print(f"[methylforge] summary -> {a.summary_out}")


def cmd_findings_xlsx(a):
    import glob
    from . import export_excel as exxl
    paths = [f for pat in a.familyforge for f in sorted(glob.glob(pat))]
    if not paths:
        raise SystemExit("no FamilyForge JSON files matched")
    exxl.cohort_from_familyforge(paths, a.out, title=a.title)
    print(f"[methylforge] cohort findings workbook ({len(paths)} samples) -> {a.out}")


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
    r.add_argument("--low-pass", action="store_true", help="relax coverage QC gates for shallow data")
    r.add_argument("--cell-reference", help="feature x cell-type reference TSV -> deconvolution")
    r.add_argument("--clock", nargs="+", help="clock coefficient TSV(s), each 'name=path' or 'path'")
    r.add_argument("--chrono-age", type=float, help="chronological age -> clock delta age")
    r.add_argument("--dmr-covariates", help="samples x covariates TSV/CSV -> adjusted DMR")
    r.add_argument("--dmr-pairs", help="sample,pair TSV/CSV -> paired DMR block")
    r.add_argument("--dmr-autosomes", action="store_true", help="drop chrX/Y/M before DMR")
    r.add_argument("--asm-hp1", help="haplotype-1 methylation file -> ASM")
    r.add_argument("--asm-hp2", help="haplotype-2 methylation file -> ASM")
    r.set_defaults(func=cmd_run)

    q = sub.add_parser("qc", help="methylation QC (conversion, coverage, sex)")
    q.add_argument("methylation")
    q.add_argument("--low-pass", action="store_true",
                   help="relax coverage gates for shallow WGBS/EM (~3-8x)")
    q.add_argument("--min-mean-cov", type=float, help="override mean-coverage gate")
    q.add_argument("--min-conversion", type=float, help="override conversion-efficiency gate")
    q.add_argument("--min-site-cov", type=int, help="override per-site min coverage")
    q.set_defaults(func=cmd_qc)

    im = sub.add_parser("imprinting"); im.add_argument("methylation")
    im.add_argument("--panel"); im.set_defaults(func=cmd_imprinting)

    d = sub.add_parser("dmr", help="differential methylation (cohort matrix), optionally covariate-adjusted/paired")
    d.add_argument("--matrix", required=True)
    d.add_argument("--groups", required=True)
    d.add_argument("--genes-bed")
    d.add_argument("--covariates", help="samples x covariates TSV/CSV (numeric or categorical) -> adjusted model")
    d.add_argument("--pair", help="sample,pair TSV/CSV -> twin/family block (paired model)")
    d.add_argument("--autosomes", action="store_true", help="drop chrX/Y/M before testing (mixed-sex contrasts)")
    d.add_argument("--model", choices=["auto", "glm"], default="auto",
                   help="'auto' two-group test, or force the covariate/paired OLS ('glm')")
    d.set_defaults(func=cmd_dmr)

    dv = sub.add_parser("dvp", help="differentially variable positions (iEVORA-style)")
    dv.add_argument("--matrix", required=True)
    dv.add_argument("--groups", required=True)
    dv.set_defaults(func=cmd_dvp)

    dc = sub.add_parser("deconvolve", help="reference-based cell-type deconvolution")
    dc.add_argument("--reference", required=True, help="feature x cell-type reference TSV (beta)")
    dc.add_argument("methylation", nargs="?", help="one per-cytosine methylation file")
    dc.add_argument("--matrix", help="feature x sample beta matrix (deconvolve every column)")
    dc.add_argument("--method", choices=["CP", "RPC"], default="CP")
    dc.add_argument("--out", help="write proportions TSV (matrix mode)")
    dc.set_defaults(func=cmd_deconvolve)

    ck = sub.add_parser("clock", help="epigenetic clock(s): DNAm age & acceleration")
    ck.add_argument("methylation")
    ck.add_argument("--clock", required=True, nargs="+",
                    help="one or more clock coefficient TSVs, each 'name=path' or 'path'")
    ck.add_argument("--chrono-age", type=float, help="chronological age -> report delta age")
    ck.set_defaults(func=cmd_clock)

    asmp = sub.add_parser("asm", help="allele-specific methylation (haplotype pair / allele table / candidate scan)")
    asmp.add_argument("--hp1", help="haplotype-1 methylation file")
    asmp.add_argument("--hp2", help="haplotype-2 methylation file")
    asmp.add_argument("--allele-table", help="long-format chrom,pos,allele,meth,cov table")
    asmp.add_argument("--candidate", help="single methylation file -> heuristic candidate-ASM scan")
    asmp.add_argument("--min-cov", type=int, default=5)
    asmp.add_argument("--min-delta", type=float, default=0.20)
    asmp.add_argument("--fdr", type=float, default=0.05)
    asmp.add_argument("--out", help="write per-CpG ASM table TSV")
    asmp.set_defaults(func=cmd_asm)

    idp = sub.add_parser("identity", help="sample-swap / identity QC from the methylation matrix")
    idp.add_argument("--matrix", help="feature x sample beta matrix")
    idp.add_argument("--samplesheet", help="CSV: sample,methylation[,group]")
    idp.add_argument("--inputs", nargs="*", help="methylation files/globs (sample = basename)")
    idp.add_argument("--expected-pairs", nargs="*", help="known-identical pairs as 'sampleA:sampleB'")
    idp.add_argument("--sex-check", action="store_true", help="add methylation-based sex per sample")
    idp.set_defaults(func=cmd_identity)

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
    co.add_argument("--autosomes-only", action="store_true", help="drop chrX/Y/M CpGs from the matrix")
    co.set_defaults(func=cmd_cohort)

    cq = sub.add_parser("cohort-qc", help="MultiQC-style batch QC dashboard across many samples")
    cq.add_argument("--samplesheet", help="CSV: sample,methylation[,group]")
    cq.add_argument("--inputs", nargs="*", help="methylation files/globs (sample = basename)")
    cq.add_argument("--results", nargs="*", help="per-sample results.json files/globs (aggregate existing QC)")
    cq.add_argument("--out", required=True, help="dashboard HTML")
    cq.add_argument("--summary-out", help="write a summary TSV")
    cq.add_argument("--title", default="MethylForge cohort QC")
    cq.set_defaults(func=cmd_cohortqc)

    fx = sub.add_parser("findings-xlsx", help="aggregate per-sample FamilyForge findings into one Excel workbook")
    fx.add_argument("--familyforge", nargs="+", required=True,
                    help="<sample>.familyforge.json files/globs")
    fx.add_argument("--out", required=True, help="output .xlsx")
    fx.add_argument("--title", default="MethylForge cohort findings")
    fx.set_defaults(func=cmd_findings_xlsx)

    t = sub.add_parser("testdata"); t.add_argument("-o", "--outdir", default="methylforge_testdata")
    t.set_defaults(func=cmd_testdata)

    v = sub.add_parser("version"); v.set_defaults(func=cmd_version)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
