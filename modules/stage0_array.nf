// Stage 0 (array) — Illumina EPIC/450K IDAT -> beta matrix + manifest (minfi/sesame).
// Downstream: `methylforge array --beta epic_beta.tsv --manifest epic_manifest.tsv`
// converts to the coordinate-keyed schema for the imprinting/cohort/episignature stages.
process ARRAY_INGEST {
  tag "epic"
  label 'array'
  publishDir "${params.outdir}/array", mode: 'copy'
  input:  path idat_dir
          path sheet
  output: path "epic_beta.tsv",     emit: beta
          path "epic_manifest.tsv", emit: manifest
          path "epic_qc.txt"
  script:
  def chain = params.array_chain ? params.array_chain : ""
  """
  Rscript ${projectDir}/bin/epic_process.R ${idat_dir} ${sheet} epic_beta.tsv epic_qc.txt ${params.array_genome} ${chain}
  """
}

// Convert the EPIC beta matrix (+manifest) into per-sample MethylForge methylation
// files, so array samples flow through the identical downstream as sequencing.
process ARRAY_TO_METHYL {
  tag "epic"
  label 'analysis'
  publishDir "${params.outdir}/array/methyl", mode: 'copy'
  input:  path beta
          path manifest
  output: path "methyl/*.methyl.tsv", emit: methyl
  script:
  """
  methylforge array --beta ${beta} --manifest ${manifest} --out-methyl-dir methyl
  """
}
