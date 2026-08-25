# Shared helper: build a probe -> GRCh38 (chrom,pos) map from the minfi 450K/EPIC
# hg19 annotations, lifted to hg38 with a UCSC chain. Self-contained (no sesame
# cache needed). Returns data.frame(probe, chrom, pos). Sourced by
# prepare_clocks.R and prepare_salas_reference.R.
build_probe_hg38 <- function(chain, platforms = c("450k", "EPIC")) {
  suppressMessages(library(rtracklayer))
  pkgs <- c()
  if ("450k" %in% platforms) pkgs <- c(pkgs, "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  if ("EPIC" %in% platforms) pkgs <- c(pkgs, "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  get_anno <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[hg38map] annotation not installed, skipping: ", pkg); return(NULL)
    }
    suppressMessages(library(pkg, character.only = TRUE))          # .onLoad assigns the object
    obj <- tryCatch(get(pkg, envir = asNamespace(pkg)), error = function(e) get(pkg))
    a <- minfi::getAnnotation(obj)
    data.frame(probe = rownames(a), chr = as.character(a$chr),
               pos = as.integer(a$pos), stringsAsFactors = FALSE)
  }
  hg19 <- do.call(rbind, lapply(pkgs, get_anno))
  if (is.null(hg19) || !nrow(hg19)) stop("no 450K/EPIC hg19 annotation available for probe mapping")
  hg19 <- hg19[!duplicated(hg19$probe) & !is.na(hg19$pos) & hg19$chr != "*", ]
  gr <- GenomicRanges::GRanges(hg19$chr, IRanges::IRanges(hg19$pos, width = 1))
  names(gr) <- hg19$probe
  lifted <- rtracklayer::liftOver(gr, rtracklayer::import.chain(chain))
  one <- lengths(lifted) == 1L
  g38 <- unlist(lifted[one])
  message(sprintf("[hg38map] %d/%d probes lifted hg19->hg38 (%s)",
                  sum(one), length(one), paste(platforms, collapse = "+")))
  data.frame(probe = names(g38), chrom = as.character(GenomicRanges::seqnames(g38)),
             pos = GenomicRanges::start(g38), stringsAsFactors = FALSE)
}

.script_dir <- function() {
  a <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("--file=", "", a[1], fixed = TRUE))) else getwd()
}
