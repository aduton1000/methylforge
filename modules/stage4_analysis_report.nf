// Stage 4-9 — MethylForge analysis + report engine (the proven Python core):
// QC, imprinting, episignature, optional DMR, HTML reports, FamilyForge export.
process ANALYSE_REPORT {
  tag "$sample"
  label 'analysis'
  publishDir "${params.outdir}/${sample}", mode: 'copy'
  input:  tuple val(sample), path(methylation)
  output: path("report/*"), emit: out
  script:
  def epi = (params.episign_ref && params.episign_labels) ? "--episign-ref ${params.episign_ref} --episign-labels ${params.episign_labels}" : ""
  def dmr = (params.run_dmr && params.dmr_matrix && params.dmr_groups) ? "--dmr-matrix ${params.dmr_matrix} --dmr-groups ${params.dmr_groups}" : ""
  def genes = params.genes_bed ? "--genes-bed ${params.genes_bed}" : ""
  """
  methylforge run ${methylation} -s ${sample} --family ${sample} \\
    --panel ${params.imprint_panel} ${epi} ${dmr} ${genes} -o report
  """
}
