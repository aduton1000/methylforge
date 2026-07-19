#!/usr/bin/env Rscript
# MethylForge — Illumina EPIC/450K IDAT -> beta matrix + probe->coordinate manifest (minfi).
# Auto-detects array type (EPIC/450K) and, for a single sample, uses noob instead of funnorm.
# For hg38 coordinates (minfi's annotation is hg19) it tries, in order: sesame's hg38
# manifest, then a deterministic UCSC liftOver of the hg19 coordinates, else keeps hg19.
# Usage:
#   Rscript epic_process.R <idat_dir> <sample_sheet.csv|auto> <out_beta.tsv> \
#           [<out_qc.txt>] [<hg38|hg19>] [<hg19ToHg38.over.chain>]
suppressMessages(library(minfi))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: epic_process.R <idat_dir> <sample_sheet|auto> <out_beta.tsv> [out_qc.txt] [hg38|hg19] [chain]")
idat_dir <- args[1]; sheet <- args[2]; out_beta <- args[3]
out_qc <- if (length(args) >= 4) args[4] else "epic_qc.txt"
genome <- if (length(args) >= 5) args[5] else "hg38"
chain  <- if (length(args) >= 6) args[6] else Sys.getenv("METHYLFORGE_CHAIN", "")

targets <- if (sheet == "auto" || !file.exists(sheet)) {
  read.metharray.sheet(idat_dir)
} else {
  read.metharray.sheet(dirname(sheet), pattern = basename(sheet))
}
rg <- read.metharray.exp(targets = targets, force = TRUE)
arraytype <- tryCatch(annotation(rg)[["array"]], error = function(e) "unknown")

detP <- detectionP(rg)
mset <- if (ncol(rg) >= 2) preprocessFunnorm(rg) else preprocessNoob(rg)
keep <- rowMeans(detP[rownames(mset), , drop = FALSE]) < 0.01
mset <- mset[keep, ]
beta <- getBeta(mset)

# hg19 coordinates from the minfi annotation
gmset <- mapToGenome(mset)
loc <- as.data.frame(SummarizedExperiment::rowRanges(gmset))
manifest <- data.frame(probe = rownames(beta), chr = as.character(loc$seqnames),
                       pos = loc$start, stringsAsFactors = FALSE)
build <- "minfi annotation (hg19)"

if (genome == "hg38") {
  remapped <- FALSE
  # 1) sesame hg38 manifest (works where sesameData's registry is current)
  if (requireNamespace("sesameData", quietly = TRUE)) {
    remapped <- tryCatch({
      plat <- if (grepl("EPIC", arraytype)) "EPIC" else "HM450"
      gr <- sesameData::sesameData_getManifestGRanges(plat, genome = "hg38")
      m2 <- as.data.frame(gr); idx <- match(manifest$probe, rownames(m2)); hit <- !is.na(idx)
      manifest <- data.frame(probe = manifest$probe[hit],
                             chr = as.character(m2$seqnames[idx[hit]]),
                             pos = m2$start[idx[hit]])
      build <- "hg38 (sesame)"; TRUE
    }, error = function(e) FALSE)
  }
  # 2) deterministic liftOver of the hg19 coordinates
  if (!remapped && nzchar(chain) && file.exists(chain) &&
      requireNamespace("rtracklayer", quietly = TRUE)) {
    remapped <- tryCatch({
      gr19 <- GenomicRanges::GRanges(manifest$chr,
                IRanges::IRanges(manifest$pos, width = 1))
      names(gr19) <- manifest$probe
      ch <- rtracklayer::import.chain(chain)
      lifted <- rtracklayer::liftOver(gr19, ch)
      one <- lengths(lifted) == 1L
      gr38 <- unlist(lifted[one])
      manifest <- data.frame(probe = names(gr38),
                             chr = as.character(GenomicRanges::seqnames(gr38)),
                             pos = GenomicRanges::start(gr38))
      build <- sprintf("hg38 (liftOver; %d/%d probes mapped)", sum(one), length(one)); TRUE
    }, error = function(e) FALSE)
  }
  if (!remapped)
    build <- "hg38 requested but no remap available -> hg19 (update sesameData or pass a chain)"
}

write.table(data.frame(probe = rownames(beta), beta, check.names = FALSE),
            file = out_beta, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(manifest, file = sub("beta", "manifest", out_beta),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("array:", arraytype), paste("genome_build:", build),
             paste("samples:", ncol(beta)), paste("probes_kept:", nrow(beta)),
             paste("probes_dropped_detP:", sum(!keep))), out_qc)
cat("wrote", out_beta, "and manifest (", build, ")\n")
