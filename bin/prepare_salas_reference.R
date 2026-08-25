#!/usr/bin/env Rscript
# MethylForge — build a cell-type deconvolution reference in MethylForge format.
#
# Produces a feature x cell-type mean-methylation TSV keyed by GRCh38 'chrom:pos'
# (the input `methylforge deconvolve --reference`) from the Salas/IDOL sorted-blood
# reference:
#   * 'extended' -> FlowSorted.BloodExtended.EPIC (Salas 2022, 12 cell types)
#   * 'blood'    -> FlowSorted.Blood.EPIC        (Salas 2018, 6 cell types)
#
# Reference *values* come from the Bioconductor data package (NOT fabricated).
# Probe -> hg38 coordinates use the sesame hg38 EPIC manifest (installed), with a
# UCSC liftOver of hg19 coordinates as a fallback (chain in _local/refs).
#
# Usage:
#   Rscript bin/prepare_salas_reference.R [extended|blood] <out.tsv> \
#           [idol|all] [<hg19ToHg38.over.chain>]
# Example:
#   Rscript bin/prepare_salas_reference.R extended \
#           _local/refs/salas_bloodextended_hg38.tsv idol _local/refs/hg19ToHg38.over.chain

suppressMessages({
  library(minfi)
})
args <- commandArgs(trailingOnly = TRUE)
.this_dir <- { a <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("--file=", "", a[1], fixed = TRUE))) else getwd() }
libname <- if (length(args) >= 1) args[1] else "extended"
out_tsv <- if (length(args) >= 2) args[2] else "salas_reference_hg38.tsv"
probeset <- if (length(args) >= 3) args[3] else "idol"
chain <- if (length(args) >= 4) args[4] else
  file.path(.this_dir, "..", "_local", "refs", "hg19ToHg38.over.chain")

pkg <- if (libname == "extended") "FlowSorted.BloodExtended.EPIC" else "FlowSorted.Blood.EPIC"
message("[salas] library = ", pkg, " ; probeset = ", probeset)
suppressMessages(library(FlowSorted.Blood.EPIC))          # provides libraryDataGet()
if (libname == "extended") suppressMessages(library(FlowSorted.BloodExtended.EPIC))

# --- reference RGChannelSet (downloaded/cached by ExperimentHub) ---------------
ref <- libraryDataGet(pkg)
message("[salas] reference loaded: ", ncol(ref), " sorted samples")

# sorted cell types only (drop whole-blood / MIX QC classes if present)
ct <- as.character(pData(ref)$CellType)
drop <- ct %in% c("MIX", "WBC", "Whole blood", "WholeBlood", "")
ref <- ref[, !drop]; ct <- ct[!drop]

# --- preprocess + per-cell-type mean beta over the IDOL probes -----------------
mset <- preprocessNoob(ref)
beta <- getBeta(mset)

idol <- NULL
if (probeset == "idol") {
  idol <- tryCatch(
    if (libname == "extended") IDOLOptimizedCpGsBloodExtended else IDOLOptimizedCpGs,
    error = function(e) NULL)
  if (is.null(idol)) message("[salas] IDOL probe list not found — using all probes")
}
probes <- if (!is.null(idol)) intersect(idol, rownames(beta)) else rownames(beta)
B <- beta[probes, , drop = FALSE]

cts <- sort(unique(ct))
ref_mat <- sapply(cts, function(c) rowMeans(B[, ct == c, drop = FALSE], na.rm = TRUE))
rownames(ref_mat) <- probes
message("[salas] reference matrix: ", nrow(ref_mat), " probes x ", ncol(ref_mat),
        " cell types (", paste(cts, collapse = ", "), ")")

# --- probe -> hg38 coordinates (450K ∪ EPIC minfi annotation + liftOver) --------
source(file.path(.this_dir, "probe_hg38_map.R"))
full <- build_probe_hg38(chain, c("450k", "EPIC"))
coords <- full[match(probes, full$probe), ]
coords <- coords[!is.na(coords$pos), ]
message("[salas] mapped ", nrow(coords), "/", length(probes), " probes to hg38")

# --- write MethylForge reference TSV (feature = chrom:pos) ----------------------
keep <- match(coords$probe, rownames(ref_mat))
out <- data.frame(feature = paste0(coords$chrom, ":", coords$pos),
                  round(ref_mat[keep, , drop = FALSE], 4), check.names = FALSE)
out <- out[!duplicated(out$feature), ]
write.table(out, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("source_package:", pkg, as.character(packageVersion(pkg))),
  paste("cell_types:", paste(cts, collapse = ",")),
  paste("probeset:", probeset, "(", length(probes), "probes )"),
  paste("features_written:", nrow(out)),
  "genome_build: GRCh38 (450K/EPIC minfi annotation + UCSC liftOver)",
  paste("built:", as.character(Sys.time()))),
  sub("\\.tsv$", ".provenance.txt", out_tsv))
cat("[salas] wrote", out_tsv, "(", nrow(out), "features x", ncol(ref_mat), "cell types )\n")
