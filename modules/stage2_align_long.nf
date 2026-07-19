// Stage 2 (long) — align ONT/PacBio native-5mC reads and pileup methylation.
// Basecalling with modified-base models (dorado, GPU) is done upstream; input is a
// modBAM (MM/ML tags) or FASTQ; output is a bedMethyl of per-CpG methylation.
process ALIGN_LONG {
  tag "$sample"
  label 'align_long'
  publishDir "${params.outdir}/${sample}/align", mode: 'copy'
  input:  tuple val(sample), path(reads)
          path genome
  output: tuple val(sample), path("${sample}.cpg.bedmethyl"), emit: methylation
  script:
  """
  # align preserving modification tags (-y), sort, index
  if samtools view -H ${reads} 2>/dev/null | grep -q '@HD'; then
     samtools fastq -T MM,ML ${reads} > reads.fq   # modBAM -> tagged fastq
     minimap2 -ax map-ont -y -t ${task.cpus} ${genome} reads.fq | samtools sort -@ ${task.cpus} -o ${sample}.aln.bam
  else
     minimap2 -ax map-ont -y -t ${task.cpus} ${genome} ${reads} | samtools sort -@ ${task.cpus} -o ${sample}.aln.bam
  fi
  samtools index ${sample}.aln.bam
  modkit pileup ${sample}.aln.bam ${sample}.cpg.bedmethyl --cpg --ref ${genome} --threads ${task.cpus}
  """
}
