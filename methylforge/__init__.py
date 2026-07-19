"""MethylForge — a reproducible epigenome (DNA-methylation) analysis pipeline.

Organism/disease-agnostic. Handles short-read (bisulfite / EM-seq) and long-read
(Nanopore / PacBio native 5mC) methylation data, genome-wide or targeted, plus
methylation-array input. Provides methylation QC, imprinting-DMR analysis,
differential-methylation (DMR/DMP) calling, episignature classification, and a
clinical/research HTML report engine.

The analysis + report engine is pure-Python (numpy/pandas/scipy/matplotlib/jinja2)
so it runs anywhere; the upstream alignment/extraction stages are Nextflow modules
that call the standard tools (Bismark/bwa-meth/MethylDackel, dorado/minimap2/modkit).
"""
__version__ = "0.1.0"
__all__ = ["io", "qc", "imprinting", "dmr", "episignature", "annotate", "plots", "report", "pipeline"]
