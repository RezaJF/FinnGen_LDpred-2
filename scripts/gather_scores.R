#!/usr/bin/env Rscript

## gather_scores.R
##
## Sum the per-chunk partial PGS files produced by score_chunk.R into one final
## per-individual PGS per weight file. Partial files are named
##   <scorebase>@@<chunk_tag>.partial.tsv   (columns: FINNGENID  partial_PGS)
## Grouping is by <scorebase> (everything before "@@"). Final PGS_i is the sum of
## partial_PGS_i over all chunks (a chunk contributes 0 where it had no matched
## variants), which reconstitutes the full genome-wide dot product.

suppressMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--partials", type = "character",
              help = "Manifest file (one partial path per line) OR a comma-separated list OR a glob directory"),
  make_option("--out_dir",  type = "character", default = "sscore",
              help = "Output directory [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

log_msg <- function(...) cat(sprintf("[gather_scores] %s\n", sprintf(...)), file = stderr())

if (is.null(opt$partials)) stop("Required: --partials")
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- resolve the list of partial files ----------------------------------
if (dir.exists(opt$partials)) {
  files <- list.files(opt$partials, pattern = "\\.partial\\.tsv(\\.gz)?$", full.names = TRUE)
} else if (file.exists(opt$partials) && !grepl(",", opt$partials)) {
  files <- readLines(opt$partials)
} else {
  files <- strsplit(opt$partials, ",")[[1]]
}
files <- trimws(files)
files <- files[nchar(files) > 0 & file.exists(files)]
if (length(files) == 0) stop("no partial-score files found")
log_msg("found %d partial files", length(files))

## ---- group by score base ------------------------------------------------
group_of <- sub("@@.*$", "", basename(files))
groups <- split(files, group_of)
log_msg("grouping into %d weight file(s)", length(groups))

## Completeness gate: every weight file must have one partial per chunk, so all
## groups must have the SAME count. A short group means partials went missing
## (dropped in flatten/localization) and the final PGS would be silently
## under-summed -- refuse rather than emit a wrong score.
grp_n <- lengths(groups)
if (length(unique(grp_n)) != 1L) {
  log_msg("uneven partial counts: %s",
          paste(sprintf("%s=%d", names(grp_n), grp_n), collapse = ", "))
  stop("gather: weight files have differing chunk counts (missing partials); ",
       "refusing to write under-summed scores")
}
log_msg("completeness OK: all %d weight file(s) have %d chunks", length(groups), grp_n[1])

qc <- vector("list", length(groups))
i <- 0
for (g in names(groups)) {
  i <- i + 1
  gfiles <- groups[[g]]

  ## Accumulate partial_PGS per FINNGENID across chunks incrementally (one FBM
  ## worth of individuals in memory, not all chunks at once). Every chunk covers
  ## the same individual set, so we keep a running total keyed by FINNGENID and
  ## align by id in case chunk row order ever differs.
  ids   <- NULL
  total <- NULL
  for (f in gfiles) {
    ## data.table in this image was compiled without zlib, so fread() cannot
    ## transparently decompress .gz; pipe through the gzip binary instead.
    src <- if (grepl("\\.gz$", f)) sprintf("gzip -dc %s", shQuote(f)) else f
    dt <- if (grepl("\\.gz$", f)) fread(cmd = src, header = TRUE, colClasses = list(character = 1L))
          else fread(f, header = TRUE, colClasses = list(character = 1L))
    setnames(dt, c("FINNGENID", "partial_PGS"))
    if (is.null(total)) {
      ids   <- dt$FINNGENID
      total <- as.numeric(dt$partial_PGS)
      names(total) <- ids
    } else if (identical(dt$FINNGENID, ids)) {
      total <- total + as.numeric(dt$partial_PGS)
    } else {
      total[dt$FINNGENID] <- total[dt$FINNGENID] + as.numeric(dt$partial_PGS)
    }
  }
  acc <- data.table(FINNGENID = ids, PGS = as.numeric(total[ids]))

  out <- file.path(opt$out_dir, paste0(g, ".sscore"))
  fwrite(acc, out, sep = "\t", quote = FALSE)
  log_msg("wrote %s (%d individuals, %d chunks summed)", out, nrow(acc), length(gfiles))

  qc[[i]] <- data.table(
    score_file      = g,
    n_individuals   = nrow(acc),
    mean            = mean(acc$PGS),
    sd              = sd(acc$PGS),
    n_chunks_summed = length(gfiles)
  )
}

qc_dt <- rbindlist(qc)
fwrite(qc_dt, file.path(opt$out_dir, "pgs_qc.tsv"), sep = "\t", quote = FALSE)
log_msg("wrote QC for %d weight file(s)", nrow(qc_dt))
