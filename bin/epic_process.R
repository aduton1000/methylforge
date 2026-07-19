#!/usr/bin/env Rscript
# MethylForge — Illumina EPIC/450K IDAT -> beta matrix + probe->coordinate manifest (minfi).
# Auto-detects array type (EPIC/450K) and, for a single sample, uses noob instead of funnorm.
# Emits hg38 coordinates via sesame when available (minfi's EPIC annotation is hg19).
# Usage: Rscript epic_process.R <idat_dir> <sample_sheet.csv|auto> <out_beta.tsv> [<out_qc.txt>] [<hg38|hg19>]
suppressMessages(library(minfi))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: epic_process.R <idat_dir> <sample_sheet|auto> <out_beta.tsv> [out_qc.txt] [hg38|hg19]")
idat_dir <- args[1]; sheet <- args[2]; out_beta <- args[3]
out_qc <- if (length(args) >= 4) args[4] else "epic_qc.txt"
genome <- if (length(args) >= 5) args[5] else "hg38"

targets <- if (sheet == "auto" || !file.exists(sheet)) {
  read.metharray.sheet(idat_dir)
} else {
  read.metharray.sheet(dirname(sheet), pattern = basename(sheet))
}
rg <- read.metharray.exp(targets = targets, force = TRUE)
arraytype <- tryCatch(annotation(rg)[["array"]], error = function(e) "unknown")

detP <- detectionP(rg)                                   # per-probe detection p-values
mset <- if (ncol(rg) >= 2) preprocessFunnorm(rg) else preprocessNoob(rg)
keep <- rowMeans(detP[rownames(mset), , drop = FALSE]) < 0.01
mset <- mset[keep, ]
beta <- getBeta(mset)

# probe -> coordinate (minfi annotation build, typically hg19)
gmset <- mapToGenome(mset)
loc <- as.data.frame(SummarizedExperiment::rowRanges(gmset))
manifest <- data.frame(probe = rownames(beta),
                       chr = as.character(loc$seqnames),
                       pos = loc$start, stringsAsFactors = FALSE)
build <- "annotation default (minfi, typically hg19)"

# best-effort hg38 remap via sesame (minfi EPIC annotation is hg19)
if (genome == "hg38" && requireNamespace("sesame", quietly = TRUE)) {
  ok <- tryCatch({
    plat <- if (grepl("EPIC", arraytype)) "EPIC" else "HM450"
    gr <- sesame::sesameData_getManifestGRanges(plat, genome = "hg38")
    m2 <- as.data.frame(gr); rn <- rownames(m2)
    idx <- match(manifest$probe, rn); hit <- !is.na(idx)
    manifest$chr[hit] <- as.character(m2$seqnames[idx[hit]])
    manifest$pos[hit] <- m2$start[idx[hit]]
    manifest <- manifest[hit, ]
    TRUE
  }, error = function(e) FALSE)
  build <- if (ok) "hg38 (sesame)" else "hg38 requested but sesame remap failed -> minfi build (hg19)"
}

write.table(data.frame(probe = rownames(beta), beta, check.names = FALSE),
            file = out_beta, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(manifest, file = sub("beta", "manifest", out_beta),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("array:", arraytype), paste("genome_build:", build),
             paste("samples:", ncol(beta)), paste("probes_kept:", nrow(beta)),
             paste("probes_dropped_detP:", sum(!keep))), out_qc)
cat("wrote", out_beta, "and manifest (", build, ")\n")
