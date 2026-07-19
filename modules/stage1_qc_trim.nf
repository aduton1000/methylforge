// Stage 1 — read QC + trimming (short-read bisulfite/EM-seq)
process QC_TRIM {
  tag "$sample"
  label 'qc_trim'
  publishDir "${params.outdir}/${sample}/qc", mode: 'copy'
  input:  tuple val(sample), path(r1), path(r2)
  output: tuple val(sample), path("${sample}_R1.trim.fq.gz"), path("${sample}_R2.trim.fq.gz"), emit: reads
          path("${sample}.fastp.json"), emit: qc
  script:
  """
  fastqc -q ${r1} ${r2} || true
  # bisulfite libraries: gentle quality/adapter trimming, keep read structure
  fastp -i ${r1} -I ${r2} -o ${sample}_R1.trim.fq.gz -O ${sample}_R2.trim.fq.gz \\
        --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 30 \\
        --json ${sample}.fastp.json --thread ${task.cpus}
  """
}
