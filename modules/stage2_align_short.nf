// Stage 2 (short) — methylation-aware alignment (Bismark default; bwa-meth alt) + dedup
process ALIGN_SHORT {
  tag "$sample"
  label 'align_short'
  publishDir "${params.outdir}/${sample}/align", mode: 'copy', pattern: "*.{bam,txt}"
  input:  tuple val(sample), path(r1), path(r2)
          path genome
  output: tuple val(sample), path("${sample}.dedup.bam"), emit: bam
  script:
  def aligner = params.aligner ?: 'bismark'
  if (aligner == 'bwameth')
  """
  bwameth.py index ${genome} 2>/dev/null || true
  bwameth.py --reference ${genome} -t ${task.cpus} ${r1} ${r2} \\
    | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam
  samtools markdup ${sample}.sorted.bam ${sample}.dedup.bam
  """
  else
  """
  GENDIR=\$(dirname ${genome})
  bismark_genome_preparation \$GENDIR 2>/dev/null || true
  bismark --genome \$GENDIR -1 ${r1} -2 ${r2} -p ${task.cpus} --basename ${sample}
  deduplicate_bismark -p --bam ${sample}_pe.bam
  mv ${sample}_pe.deduplicated.bam ${sample}.dedup.bam
  """
}
