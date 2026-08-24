#!/usr/bin/env Rscript
# MethylForge — build epigenetic-clock coefficient files in MethylForge format.
#
# Converts authoritative clock coefficients from the `dnaMethyAge` package into
# per-clock TSVs (`feature`, `coefficient`, optional `mean`) + `.meta.json`
# transform specs, and writes/refreshes resources/clocks/manifest.tsv. Coefficient
# *values* come from dnaMethyAge (NOT fabricated). Two keyspaces are emitted:
#   * <clock>.probe.tsv  — keyed by EPIC/450K probe id (for array betas, or WGBS
#                          after `methylforge array`)
#   * <clock>.hg38.tsv   — keyed by GRCh38 'chrom:pos' (for sequencing), via the
#                          sesame HM450/EPIC hg38 manifests.
#
# Usage:
#   Rscript bin/prepare_clocks.R [--out resources/clocks] [--clocks id,id,...] [--no-coords]

suppressMessages({ library(dnaMethyAge) })
args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
outdir <- getopt("--out", "_local/refs/clocks")
do_coords <- !("--no-coords" %in% args)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1], fixed = TRUE)),
                 "probe_hg38_map.R"))
chain <- getopt("--chain", file.path(.script_dir(), "..", "_local", "refs", "hg19ToHg38.over.chain"))

# clock registry: dnaMethyAge dataset id -> MethylForge name + transform metadata
# (facts from the source publications; only Horvath-family clocks use the anti-log
#  age transform with the adult_age=20 knot). DunedinPACE is intentionally excluded
# here — it needs dnaMethyAge's own quantile normalisation, not an additive clock.
REG <- list(
  HorvathS2013 = list(name="Horvath2013",          transform="horvath", adult_age=20, unit="years", tissue="multi-tissue",           pediatric=FALSE, doi="10.1186/gb-2013-14-10-r115"),
  HorvathS2018 = list(name="HorvathSkinBlood2018",  transform="horvath", adult_age=20, unit="years", tissue="skin+blood",            pediatric=FALSE, doi="10.18632/aging.101508"),
  HannumG2013  = list(name="Hannum2013",            transform="none",    adult_age=NA, unit="years", tissue="blood",                 pediatric=FALSE, doi="10.1016/j.molcel.2012.10.016"),
  LevineM2018  = list(name="PhenoAge2018",          transform="none",    adult_age=NA, unit="years", tissue="blood",                 pediatric=FALSE, doi="10.18632/aging.101414"),
  McEwenL2019  = list(name="PedBE2019",             transform="horvath", adult_age=20, unit="years", tissue="buccal (pediatric)",    pediatric=TRUE,  doi="10.1073/pnas.1820843116"),
  ZhangQ2019   = list(name="Zhang2019",             transform="none",    adult_age=NA, unit="years", tissue="blood/saliva",          pediatric=FALSE, doi="10.1186/s13073-019-0667-1")
)
want <- getopt("--clocks", "")
ids <- if (nzchar(want)) strsplit(want, ",")[[1]] else names(REG)

# probe -> hg38 coordinate map (450K ∪ EPIC minfi annotation + liftOver)
coord_map <- NULL
if (do_coords) {
  coord_map <- tryCatch(build_probe_hg38(chain, c("450k", "EPIC")),
    error = function(e) { message("[clocks] coord map unavailable (", conditionMessage(e),
                                  ") — writing probe-keyed clocks only"); NULL })
  if (!is.null(coord_map)) message("[clocks] coord map: ", nrow(coord_map), " probes -> hg38")
}

# Horvath gold-standard means (for imputing missing clock CpGs), if available
gsm <- tryCatch({ e <- new.env(); data("probeAnnotation21kdatMethUsed", package="dnaMethyAge", envir=e)
  a <- get(ls(e)[1], envir=e); mc <- grep("gold", colnames(a), ignore.case=TRUE, value=TRUE)[1]
  nc <- grep("name|cpg|probe", colnames(a), ignore.case=TRUE, value=TRUE)[1]
  if (!is.na(mc) && !is.na(nc)) setNames(a[[mc]], a[[nc]]) else NULL }, error=function(e) NULL)

manifest <- list()
for (id in ids) {
  meta <- REG[[id]]; if (is.null(meta)) { message("[clocks] skip unknown id: ", id); next }
  e <- new.env()
  ok <- tryCatch({ data(list=id, package="dnaMethyAge", envir=e); TRUE }, error=function(x) FALSE)
  if (!ok || !exists("coefs", envir=e)) { message("[clocks] not found in dnaMethyAge: ", id); next }
  co <- get("coefs", envir=e)
  colnames(co)[1:2] <- c("feature", "coefficient")
  is_int <- grepl("intercept", co$feature, ignore.case=TRUE)
  n_cpg <- sum(!is_int)
  df <- co[, c("feature", "coefficient")]
  if (!is.null(gsm) && meta$name == "Horvath2013")
    df$mean <- ifelse(is_int, NA, unname(gsm[df$feature]))

  base <- file.path(outdir, meta$name)
  write.table(df, paste0(base, ".probe.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  metaj <- sprintf('{"name": "%s", "transform": "%s", "adult_age": %s, "unit": "%s", "tissue": "%s", "source": "dnaMethyAge:%s", "doi": "%s", "keyspace": "probe"}',
                   meta$name, meta$transform, ifelse(is.na(meta$adult_age), "null", meta$adult_age),
                   meta$unit, meta$tissue, id, meta$doi)
  writeLines(metaj, paste0(base, ".probe.tsv.meta.json"))
  n_hg38 <- NA
  if (!is.null(coord_map)) {
    idx <- match(df$feature, coord_map$probe)
    dfc <- df; dfc$feature <- ifelse(is_int, "(Intercept)",
              ifelse(is.na(idx), NA, paste0(coord_map$chrom[idx], ":", coord_map$pos[idx])))
    dfc <- dfc[is_int | !is.na(dfc$feature), ]
    n_hg38 <- sum(!grepl("Intercept", dfc$feature))
    write.table(dfc, paste0(base, ".hg38.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    writeLines(sub('"keyspace": "probe"', '"keyspace": "chrom:pos"', metaj),
               paste0(base, ".hg38.tsv.meta.json"))
  }
  manifest[[length(manifest)+1]] <- data.frame(
    clock=meta$name, source_id=id, transform=meta$transform,
    adult_age=ifelse(is.na(meta$adult_age), "", meta$adult_age),
    n_cpg=n_cpg, n_cpg_hg38=ifelse(is.na(n_hg38), "", n_hg38),
    unit=meta$unit, tissue=meta$tissue, pediatric=meta$pediatric, doi=meta$doi,
    stringsAsFactors=FALSE)
  message(sprintf("[clocks] %-22s %d CpGs (%s hg38-mapped) transform=%s",
                  meta$name, n_cpg, ifelse(is.na(n_hg38), "-", n_hg38), meta$transform))
}
if (length(manifest)) {
  man <- do.call(rbind, manifest)
  write.table(man, file.path(outdir, "manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  cat("[clocks] wrote", nrow(man), "clocks +", file.path(outdir, "manifest.tsv"), "\n")
  cat("[clocks] NOTE: DunedinPACE excluded (needs dnaMethyAge::preprocessDunedinPACE normalisation).\n")
} else cat("[clocks] no clocks built\n")
