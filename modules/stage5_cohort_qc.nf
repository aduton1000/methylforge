// Cohort QC dashboard — aggregate per-sample methylation QC across the run into
// one MultiQC-style HTML + summary TSV.
process COHORT_QC {
  tag "cohort"
  label 'report'
  publishDir "${params.outdir}", mode: 'copy'
  input:  path methfiles
  output: path "cohort_qc.html"
          path "cohort_qc_summary.tsv"
  script:
  """
  methylforge cohort-qc --inputs ${methfiles} \\
    --out cohort_qc.html --summary-out cohort_qc_summary.tsv \\
    --title "MethylForge cohort QC"
  """
}
