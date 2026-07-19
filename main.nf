#!/usr/bin/env nextflow
// MethylForge — main entry. Select a mode with --mode:
//   --mode sequencing  (default) short/long/precalled -> analysis + reports + batch QC
//   --mode array       Illumina EPIC/450K IDAT -> beta (minfi) -> analysis + reports
//   --mode cohort_qc   batch QC dashboard only, from a samplesheet of methylation files
// e.g. nextflow run main.nf -profile local,conda --mode array --idat_dir DIR --array_sheet s.csv
nextflow.enable.dsl = 2

include { METHYLFORGE } from './subworkflows/methylforge.nf'
include { ARRAY_INGEST; ARRAY_TO_METHYL } from './modules/stage0_array.nf'
include { COHORT_QC } from './modules/stage5_cohort_qc.nf'

def parse_row(row) {
  def meta = [ id: row.sample,
               read_type: (row.read_type ?: params.read_type),
               group: (row.group ?: 'NA') ]
  def files = [:]
  if (meta.read_type == 'short')     { files.r1 = file(row.fastq_1); files.r2 = file(row.fastq_2) }
  else if (meta.read_type == 'long') { files.reads = file(row.modbam ?: row.fastq_1) }
  else                               { files.methylation = file(row.methylation) }
  return tuple(meta, files)
}

// ── sequencing / precalled ────────────────────────────────────────────────
workflow sequencing {
  if (!params.samplesheet) error "provide --samplesheet (sample,fastq_1,fastq_2,read_type,modbam,methylation,group)"
  samples = Channel.fromPath(params.samplesheet).splitCsv(header: true).map { parse_row(it) }
  genome  = params.reference ? file(params.reference) : file("${projectDir}/assets/NO_FILE")
  METHYLFORGE(samples, genome)
}

// ── array: IDAT -> beta (minfi) -> per-sample methylation -> analysis ──────
//   from IDATs:            --idat_dir DIR --array_sheet sheet.csv
//   from a pre-made beta:  --array_beta beta.tsv --array_manifest manifest.tsv (skips the R step)
workflow array {
  if (params.array_beta && params.array_manifest) {
    ARRAY_TO_METHYL(file(params.array_beta), file(params.array_manifest))
  } else {
    if (!params.idat_dir || !params.array_sheet)
      error "array mode needs --idat_dir + --array_sheet, or --array_beta + --array_manifest"
    ARRAY_INGEST(file(params.idat_dir), file(params.array_sheet))
    ARRAY_TO_METHYL(ARRAY_INGEST.out.beta, ARRAY_INGEST.out.manifest)
  }
  samples = ARRAY_TO_METHYL.out.methyl.flatten()
      .map { f -> tuple([id: f.getBaseName().replaceAll(/\.methyl$/, ''),
                         read_type: 'precalled', group: 'NA'], [methylation: f]) }
  METHYLFORGE(samples, file("${projectDir}/assets/NO_FILE"))
}

// ── cohort_qc: batch QC dashboard only ────────────────────────────────────
workflow cohort_qc {
  if (!params.samplesheet) error "cohort_qc mode needs --samplesheet (sample,methylation)"
  meth = Channel.fromPath(params.samplesheet).splitCsv(header: true)
                .map { file(it.methylation) }.collect()
  COHORT_QC(meth)
}

// ── dispatch ──────────────────────────────────────────────────────────────
workflow {
  if      (params.mode == 'array')     { array() }
  else if (params.mode == 'cohort_qc') { cohort_qc() }
  else                                 { sequencing() }
}
