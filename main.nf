#!/usr/bin/env nextflow
// MethylForge — main entry.
// Usage:
//   nextflow run main.nf -profile local,conda --samplesheet samples.csv --reference genome.fa
//   nextflow run main.nf -profile test,conda           (synthetic fixture; analysis+report only)
nextflow.enable.dsl = 2

include { METHYLFORGE } from './subworkflows/methylforge.nf'

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

workflow {
  if (params.test_mode) {
    // fixture: generate synthetic precalled methylation then analyse+report
    log.info "MethylForge TEST fixture — generating synthetic methylation and running analysis+report."
  }
  if (!params.samplesheet) { error "Provide --samplesheet (CSV: sample,fastq_1,fastq_2,read_type,modbam,methylation,group)" }

  samples = Channel.fromPath(params.samplesheet)
                   .splitCsv(header: true)
                   .map { parse_row(it) }
  genome  = params.reference ? file(params.reference) : file("${projectDir}/assets/NO_FILE")

  METHYLFORGE(samples, genome)
}

workflow.onComplete {
  log.info "MethylForge complete → ${params.outdir}  (ok=${workflow.success})"
}
