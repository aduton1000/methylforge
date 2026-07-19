// MethylForge subworkflow — route each sample from its input type to a normalised
// per-CpG methylation file, then run the analysis + report engine.
include { QC_TRIM }        from '../modules/stage1_qc_trim.nf'
include { ALIGN_SHORT }    from '../modules/stage2_align_short.nf'
include { ALIGN_LONG }     from '../modules/stage2_align_long.nf'
include { EXTRACT }        from '../modules/stage3_extract.nf'
include { ANALYSE_REPORT } from '../modules/stage4_analysis_report.nf'
include { COHORT_QC }      from '../modules/stage5_cohort_qc.nf'

workflow METHYLFORGE {
  take:
    samples          // channel: tuple(meta, files...) with meta.read_type in short|long|precalled
    genome

  main:
    // split by read type
    short_ch  = samples.filter { it[0].read_type == 'short'   }.map { m, f -> tuple(m.id, f.r1, f.r2) }
    long_ch   = samples.filter { it[0].read_type == 'long'    }.map { m, f -> tuple(m.id, f.reads) }
    pre_ch    = samples.filter { it[0].read_type == 'precalled' }.map { m, f -> tuple(m.id, f.methylation) }

    // short-read path: QC/trim -> align -> extract
    QC_TRIM(short_ch)
    ALIGN_SHORT(QC_TRIM.out.reads, genome)
    EXTRACT(ALIGN_SHORT.out.bam, genome)

    // long-read path: align + pileup
    ALIGN_LONG(long_ch, genome)

    // converge all methylation call sets
    methylation = EXTRACT.out.methylation
                    .mix(ALIGN_LONG.out.methylation)
                    .mix(pre_ch)

    ANALYSE_REPORT(methylation)

    // batch QC dashboard across all samples in the run
    if (params.cohort_qc)
        COHORT_QC(methylation.map { s, f -> f }.collect())

  emit:
    reports = ANALYSE_REPORT.out.out
}
