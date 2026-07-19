#!/usr/bin/env Rscript
# MethylForge — Illumina EPIC/450K IDAT -> beta matrix + manifest (minfi).
# Usage: Rscript epic_process.R <idat_dir> <sample_sheet.csv> <out_beta.tsv> [<out_qc.txt>]
suppressMessages(library(minfi))
args <- commandArgs(trailingOnly = TRUE)
idat_dir <- args[1]; sheet <- args[2]; out_beta <- args[3]
out_qc <- ifelse(length(args) >= 4, args[4], "epic_qc.txt")

targets <- read.metharray.sheet(idat_dir, pattern = sheet)
rg <- read.metharray.exp(targets = targets, force = TRUE)
detP <- detectionP(rg)                       # per-probe detection p-values (QC)
mset <- preprocessFunnorm(rg)                 # functional normalisation (use preprocessNoob for small sets)
keep <- rowMeans(detP[rownames(mset), , drop = FALSE]) < 0.01
mset <- mset[keep, ]
beta <- getBeta(mset)

write.table(data.frame(probe = rownames(beta), beta, check.names = FALSE),
            file = out_beta, sep = "\t", quote = FALSE, row.names = FALSE)
ann <- getAnnotation(mset)[, c("chr", "pos")]  # probe -> genomic coordinate
write.table(data.frame(probe = rownames(ann), chr = ann$chr, pos = ann$pos),
            file = sub("beta", "manifest", out_beta), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("samples:", ncol(beta)), paste("probes_kept:", nrow(beta)),
             paste("probes_dropped_detP:", sum(!keep))), out_qc)
cat("wrote", out_beta, "and manifest\n")
