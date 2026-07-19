// Stage 3 — methylation extraction from short-read BAM -> per-CpG methylation.
process EXTRACT {
  tag "$sample"
  label 'extract'
  publishDir "${params.outdir}/${sample}/methylation", mode: 'copy'
  input:  tuple val(sample), path(bam)
          path genome
  output: tuple val(sample), path("${sample}.cpg.bedGraph"), emit: methylation
  script:
  """
  samtools sort -@ ${task.cpus} -o ${sample}.srt.bam ${bam}
  samtools index ${sample}.srt.bam
  samtools faidx ${genome}
  # per-CpG methylation; --CHH/--CHG enable conversion-rate QC downstream
  MethylDackel extract --CHH --CHG -@ ${task.cpus} ${genome} ${sample}.srt.bam -o ${sample}
  # MethylDackel writes ${sample}_CpG.bedGraph
  mv ${sample}_CpG.bedGraph ${sample}.cpg.bedGraph
  """
}
