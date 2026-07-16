#!/usr/bin/env Rscript

## combine_ld_reference.R  (optional convenience)
##
## Bundles the per-chromosome LD parts from build_ld_reference.R into a SINGLE
## portable .rds, so a precomputed LD reference can be distributed as one file:
##     list(parts = list(<part1>, <part2>, ...))       # chr-ordered
## where each part is list(corr = <dsCMatrix>, map = data.frame(chr,pos,a0,a1),
## ld = <numeric>).
##
## This is NOT required by the pipeline: ldpred2_infer.R can consume the per-chr
## part files directly (--ld_parts manifest). Use this only to package a
## reference for reuse/sharing. The bundle stays fully portable (plain sparse
## matrices; the on-disk SFBM is (re)built inside ldpred2_infer.R at run time).

suppressMessages({ library(optparse); library(data.table) })
log_msg <- function(...) cat(sprintf("[combine_ld] %s\n", sprintf(...)), file = stderr())

opt <- parse_args(OptionParser(option_list = list(
  make_option("--parts", type = "character", help = "manifest (one per-chr .rds path/line) or comma list"),
  make_option("--out",   type = "character", help = "output single-file LD bundle .rds")
)))
for (r in c("parts","out")) if (is.null(opt[[r]])) stop("Required: --parts --out")

if (file.exists(opt$parts) && !grepl(",", opt$parts)) paths <- readLines(opt$parts) else paths <- strsplit(opt$parts, ",")[[1]]
paths <- trimws(paths); paths <- paths[nchar(paths) > 0 & file.exists(paths)]
if (length(paths) == 0) stop("no per-chr LD parts found")

parts <- lapply(paths, readRDS)
for (p in parts) if (is.null(p$corr) || is.null(p$map) || is.null(p$ld))
  stop("each part must be list(corr, map, ld)")
## order by chromosome for a deterministic bundle (X/Y/MT -> 23..26, like the rest of the pipeline)
chr_to_int <- function(x) {
  x <- sub("^CHR", "", toupper(as.character(x)))
  x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x == "XY"] <- "25"; x[x %in% c("MT","M")] <- "26"
  suppressWarnings(as.integer(x))
}
chr1 <- sapply(parts, function(p) { m <- as.data.frame(p$map); if (nrow(m)) chr_to_int(m$chr[1]) else NA_integer_ })
parts <- parts[order(chr1)]

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(parts = parts), opt$out)
log_msg("bundled %d parts (%d variants) -> %s",
        length(parts), sum(sapply(parts, function(p) nrow(as.data.frame(p$map)))), opt$out)
