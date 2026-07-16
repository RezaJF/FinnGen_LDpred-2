#!/usr/bin/env Rscript

## score_chunk.R
##
## Pure PGS SCORING for ONE FinnGen R13 bgen chunk using bigsnpr.
## Applies precomputed (already inferred) per-SNP weights to imputed dosages.
## This is a genotype x weight dot product: PGS_i = sum_j dosage_ij * beta_j.
## There is NO LDpred2 inference and NO LD matrix used here.
##
## For each harmonised weight file it emits one partial-score file:
##   <scorebase>@@<chunk_tag>.partial.tsv  with columns  FINNGENID  partial_PGS
## gather_scores.R then sums these partials across all chunks per weight file.
##
## bigsnpr functions used (argument names per bigsnpr >= 1.9, verify against the
## pinned image version):
##   snp_readBGEN(bgenfiles, backingfile, list_snp_id, ncores)
##       - list_snp_id is a LIST with one character vector per bgen file, each id
##         in FinnGen bgen format "<chr>_<pos>_<a1>_<a2>". EVERY requested id must
##         exist in the bgen or snp_readBGEN errors, so we first intersect the
##         requested ids with the ids actually present in this chunk (read from
##         the .bgi SQLite index) before calling it.
##   snp_attach(rds)                  - load the produced bigSNP object.
##   snp_match(sumstats, info_snp)    - align weights to chunk map, flipping the
##                                      beta sign for allele swaps / reversals.
##   big_prodVec(G, y.col, ind.col)   - the actual G[, ind] %*% beta dot product.
##
## Dosage/allele convention (verify on first run): snp_readBGEN returns the
## dosage of the SECOND allele of each "chr_pos_a1_a2" id, and obj$map stores
## allele1 (first) / allele2 (second). We build info_snp with a0 = allele1 and
## a1 = allele2, so snp_match returns beta with respect to a1 = the counted
## dosage allele, making big_prodVec directly correct.

suppressMessages({
  library(optparse)
  library(data.table)
  library(bigsnpr)
  library(bigreadr)
  library(RSQLite)
})

## bigstatsr aborts with "Two levels of parallelism are used" when its own
## ncores>1 parallelism runs on top of a multi-threaded BLAS. Pin BLAS to a
## single thread (this targets ONLY the BLAS, leaving bigsnpr/bigstatsr's ncores
## compute parallelism intact) so big_prodVec/snp_readBGEN can use --ncores.
suppressWarnings(try(bigparallelr::set_blas_ncores(1), silent = TRUE))

option_list <- list(
  make_option("--bgen",         type = "character", help = "Path to one .bgen chunk"),
  make_option("--sample",       type = "character", help = "Path to the .bgen.sample (Oxford) file"),
  make_option("--score_files",  type = "character", help = "Comma-separated harmonised score files, or a manifest file (one path per line)"),
  make_option("--id_col",       type = "character", default = "bgen_id",   help = "Column holding the FinnGen bgen variant id [default %default]"),
  make_option("--a1_col",       type = "character", default = "a1_effect", help = "Effect-allele column [default %default]"),
  make_option("--a0_col",       type = "character", default = "a0",        help = "Other-allele column [default %default]"),
  make_option("--weight_col",   type = "character", default = "weight",    help = "Per-SNP weight/beta column [default %default]"),
  make_option("--chunk_tag",    type = "character", default = NULL,        help = "Tag to make partial filenames unique [default: bgen basename]"),
  make_option("--ncores",       type = "integer",   default = 1,           help = "Threads for bigsnpr [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

log_msg <- function(...) cat(sprintf("[score_chunk] %s\n", sprintf(...)), file = stderr())

if (is.null(opt$bgen) || is.null(opt$sample) || is.null(opt$score_files)) {
  stop("Required: --bgen, --sample, --score_files")
}
chunk_tag <- if (is.null(opt$chunk_tag)) sub("\\.bgen$", "", basename(opt$bgen)) else opt$chunk_tag

## ---- resolve the list of score files ------------------------------------
if (file.exists(opt$score_files) && !grepl(",", opt$score_files)) {
  score_paths <- readLines(opt$score_files)
} else {
  score_paths <- strsplit(opt$score_files, ",")[[1]]
}
score_paths <- trimws(score_paths)
score_paths <- score_paths[nchar(score_paths) > 0]
log_msg("scoring %d weight file(s) against chunk %s", length(score_paths), chunk_tag)

scorebase <- function(p) {
  b <- basename(p)
  b <- sub("\\.gz$", "", b)
  b <- sub("\\.score\\.tsv$", "", b)
  b <- sub("\\.tsv$", "", b)
  b
}

## Partial-score filenames (and gather's grouping) key on scorebase(), so two
## weight files reducing to the same basename would silently clobber each other.
bases <- vapply(score_paths, scorebase, character(1))
if (anyDuplicated(bases)) {
  dup <- unique(bases[duplicated(bases)])
  stop(sprintf("duplicate score-file basename(s) would collide: %s",
               paste(dup, collapse = ", ")))
}

## ---- helper: normalise chromosome to integer (1..22, X=23, Y=24, MT=26) --
chr_to_int <- function(x) {
  x <- sub("^CHR", "", toupper(as.character(x)))
  x[x == "X"]  <- "23"
  x[x == "Y"]  <- "24"
  x[x == "XY"] <- "25"
  x[x %in% c("MT", "M")] <- "26"
  suppressWarnings(as.integer(x))
}

## ---- read the harmonised weight files -----------------------------------
score_dt <- list()
all_ids  <- character(0)
for (p in score_paths) {
  dt <- fread(p)
  need <- c(opt$id_col, opt$a1_col, opt$a0_col, opt$weight_col)
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop(sprintf("score file %s missing columns: %s", p, paste(miss, collapse = ", ")))
  dt <- dt[!is.na(get(opt$weight_col)) & get(opt$weight_col) != 0]
  score_dt[[p]] <- dt
  all_ids <- union(all_ids, as.character(dt[[opt$id_col]]))
}
log_msg("union of requested variant ids across weight files: %d", length(all_ids))

## ---- which of those variants are actually in THIS chunk ------------------
## snp_readBGEN errors if a requested id is absent, so intersect first using the
## .bgi SQLite index. We reconstruct ids as "<chr>_<pos>_<allele1>_<allele2>"
## from the index so the format matches what snp_readBGEN expects internally.
bgi_path <- paste0(opt$bgen, ".bgi")
con <- dbConnect(RSQLite::SQLite(), bgi_path)
on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)
vt <- as.data.table(dbGetQuery(con, "SELECT chromosome, position, allele1, allele2 FROM Variant"))
try(dbDisconnect(con), silent = TRUE)

## CRITICAL: bigsnpr::snp_readBGEN matches list_snp_id via
##   match(format_snp_id(list_snp_id), format_snp_id(<chr>_<pos>_<a1>_<a2> from the .bgi))
## and format_snp_id() ONLY accepts a numeric chromosome ("1" -> "01"); it aborts
## with "Wrong format of some variants" on FinnGen's "chr1". So we rewrite the
## localized .bgi 'chromosome' column to a bare integer IN PLACE (the Cromwell-
## localized copy is writable and isolated per shard), and build the requested ids
## with the same bare-integer chromosome. snp_readBGEN then reads that rewritten
## index (info_id numeric) and our numeric list_snp_id matches it. Only the text
## column is touched; file_start_position/size_in_bytes (the seek offsets) are
## untouched, so dosage reads are unaffected.
chrmap <- unique(vt$chromosome)
con_up <- dbConnect(RSQLite::SQLite(), bgi_path)
for (ch in chrmap) {
  dbExecute(con_up, "UPDATE Variant SET chromosome = :new WHERE chromosome = :old",
            params = list(new = as.character(chr_to_int(ch)), old = ch))
}
dbDisconnect(con_up)
## numeric-chromosome ids, matching the rewritten .bgi that snp_readBGEN will read
bgi_ids <- paste(chr_to_int(vt$chromosome), vt$position, vt$allele1, vt$allele2, sep = "_")

## Match harmonised weights to chunk variants on a CANONICAL key
##   chr-as-int : pos : {allele pair, order-independent}
## so the "which variants to read" decision is robust to allele ORDER
## (ref/alt vs allele1/allele2) and to chromosome-string convention ("chr1" vs
## "1"). snp_match() below still reconciles the effect-allele orientation and
## flips beta; here we only pick which chunk variants to READ. want_ids below
## uses the bare-integer .bgi id (built above) that the rewritten index and
## snp_readBGEN's format_snp_id() both accept.
allele_key <- function(a, b) {
  a <- toupper(as.character(a)); b <- toupper(as.character(b))
  ifelse(a < b, paste(a, b, sep = "/"), paste(b, a, sep = "/"))
}
bgi_key <- paste(chr_to_int(vt$chromosome), vt$position,
                 allele_key(vt$allele1, vt$allele2), sep = ":")

harm_keys <- character(0)
for (dt in score_dt) {
  if (all(c("chr", "pos") %in% names(dt))) {
    k <- paste(chr_to_int(dt[["chr"]]), dt[["pos"]],
               allele_key(dt[[opt$a1_col]], dt[[opt$a0_col]]), sep = ":")
  } else {
    parts <- tstrsplit(as.character(dt[[opt$id_col]]), "_", fixed = TRUE)
    k <- paste(chr_to_int(parts[[1]]), parts[[2]],
               allele_key(parts[[3]], parts[[4]]), sep = ":")
  }
  harm_keys <- union(harm_keys, k)
}

want_ids <- unique(bgi_ids[bgi_key %in% harm_keys])
log_msg("variants requested present in this chunk (allele/chr-robust match): %d of %d union",
        length(want_ids), length(all_ids))
if (length(want_ids) == 0L)
  log_msg("WARNING: zero variants matched in this chunk - check chr convention / allele encoding")

n_ind <- NA_integer_
finngenid <- NULL

## ---- read the .sample file to recover FINNGENID order --------------------
## Oxford .sample: line 1 header, line 2 types, then one row per individual in
## bgen genotype-row order. FinnGen files carry the id in column ID_2 (col 2);
## fall back to column 1 if column 2 is all "0". VERIFY column choice on first run.
samp <- fread(opt$sample, header = FALSE, skip = 2)
if (ncol(samp) >= 2 && !all(as.character(samp[[2]]) == "0")) {
  finngenid <- as.character(samp[[2]])
} else {
  finngenid <- as.character(samp[[1]])
}
n_ind <- length(finngenid)
log_msg("individuals in chunk (from .sample): %d", n_ind)

write_partial <- function(base, values) {
  # Gzipped to keep gather-stage localization tractable across 32 x 535
  # partials of ~520k rows each. data.table in this image was compiled
  # without zlib, so fwrite() cannot gzip by extension; write plain then
  # gzip via the system binary (rocker base ships gzip).
  out   <- sprintf("%s@@%s.partial.tsv.gz", base, chunk_tag)
  plain <- sub("\\.gz$", "", out)
  fwrite(data.table(FINNGENID = finngenid, partial_PGS = values),
         plain, sep = "\t", quote = FALSE)
  rc <- system2("gzip", c("-f", shQuote(plain)))
  if (rc != 0L || !file.exists(out))
    stop(sprintf("gzip failed for %s (rc=%s)", plain, rc))
  log_msg("wrote %s", out)
}

## ---- defensive: no overlapping variants in this chunk at all -------------
if (length(want_ids) == 0) {
  log_msg("WARNING: zero requested variants present in chunk %s; writing zeros", chunk_tag)
  for (p in score_paths) write_partial(scorebase(p), rep(0, n_ind))
  quit(save = "no", status = 0)
}

## ---- read the needed dosages once for the whole chunk --------------------
backing <- tempfile(tmpdir = getwd(), fileext = "")
rds <- snp_readBGEN(
  bgenfiles    = opt$bgen,
  backingfile  = backing,
  list_snp_id  = list(want_ids),   # LIST: one character vector per bgen file
  ncores       = opt$ncores
)
obj <- snp_attach(rds)
G   <- obj$genotypes
map <- obj$map
# map columns: chromosome, marker.ID, genetic.dist, physical.pos, allele1, allele2
log_msg("bigSNP loaded: %d individuals x %d variants", nrow(G), ncol(G))

if (nrow(G) != n_ind) {
  log_msg("WARNING: bgen individual count (%d) != .sample rows (%d); trusting bgen order",
          nrow(G), n_ind)
}

## info_snp for snp_match: a0 = allele1 (first), a1 = allele2 (counted dosage).
info_snp <- data.frame(
  chr = chr_to_int(map$chromosome),
  pos = as.integer(map$physical.pos),
  a0  = as.character(map$allele1),
  a1  = as.character(map$allele2),
  stringsAsFactors = FALSE
)

## ---- score each weight file ---------------------------------------------
for (p in score_paths) {
  base <- scorebase(p)
  dt <- score_dt[[p]]

  sumstats <- data.frame(
    chr  = chr_to_int(sub("_.*$", "", as.character(dt[[opt$id_col]]))),
    pos  = as.integer(sub("^[^_]*_([^_]*)_.*$", "\\1", as.character(dt[[opt$id_col]]))),
    a0   = as.character(dt[[opt$a0_col]]),
    a1   = as.character(dt[[opt$a1_col]]),
    beta = as.numeric(dt[[opt$weight_col]]),
    stringsAsFactors = FALSE
  )
  sumstats <- sumstats[!is.na(sumstats$chr) & !is.na(sumstats$pos), , drop = FALSE]

  # strand_flip = FALSE: the weights are already harmonised to the target bgen's
  # build and alleles (same cohort), so no strand-complement resolution is needed.
  # Leaving the default (TRUE) would DROP every strand-ambiguous variant (A/T, C/G)
  # -- ~1/3 of SNPs -- silently biasing each score. Reversal (a0/a1 swap) is still
  # handled independently of this flag.
  matched <- tryCatch(
    snp_match(sumstats, info_snp, join_by_pos = TRUE, match.min.prop = 0,
              strand_flip = FALSE),
    error = function(e) { log_msg("snp_match failed for %s: %s", base, conditionMessage(e)); NULL }
  )

  if (is.null(matched) || nrow(matched) == 0) {
    log_msg("%s: 0 variants matched in chunk %s -> writing zeros", base, chunk_tag)
    write_partial(base, rep(0, nrow(G)))
    next
  }

  n_flip <- sum(matched$beta != sumstats$beta[match(paste(matched$chr, matched$pos),
                                                    paste(sumstats$chr, sumstats$pos))], na.rm = TRUE)
  log_msg("%s: matched %d / %d weights (chunk %s); ~%d sign-flipped",
          base, nrow(matched), nrow(sumstats), chunk_tag, n_flip)

  ind_col <- matched[["_NUM_ID_"]]   # column index into G
  beta    <- matched$beta            # aligned to a1 = counted dosage allele

  partial <- big_prodVec(G, beta, ind.col = ind_col, ncores = opt$ncores)
  if (anyNA(partial)) {
    log_msg("%s: %d NA partial values -> set to 0 (unexpected for dosages)",
            base, sum(is.na(partial)))
    partial[is.na(partial)] <- 0
  }
  write_partial(base, partial)
}

log_msg("chunk %s done", chunk_tag)
