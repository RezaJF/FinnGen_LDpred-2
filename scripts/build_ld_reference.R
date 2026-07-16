#!/usr/bin/env Rscript

## build_ld_reference.R  (per-chromosome)
##
## Computes an in-sample LD (correlation) reference for ONE chromosome from a
## bgen chunk, on a subsample of individuals, restricted to a supplied variant
## list (e.g. HapMap3+). Emits a per-chromosome bundle (sparse corr + map + LD
## scores) that combine_ld_reference.R assembles into the genome-wide reference
## consumed by ldpred2_infer.R.
##
## LD is estimated on a random subsample of individuals (default 10,000) because
## the correlation reference does not need the full cohort and this keeps memory
## and time bounded (standard practice for LDpred2).
##
## Uses the same FinnGen bgen conventions as score_chunk.R: reads the .bgi
## SQLite index, rewrites its 'chromosome' text column from 'chr1' to a bare
## integer in place (bigsnpr::format_snp_id needs numeric chromosomes), and
## matches variants on a canonical chr:pos:{sorted alleles} key.

suppressMessages({
  library(optparse)
  library(data.table)
  library(bigsnpr)
  library(RSQLite)
  library(Matrix)
})
log_msg <- function(...) cat(sprintf("[build_ld] %s\n", sprintf(...)), file = stderr())

option_list <- list(
  make_option("--bgen",        type = "character", help = "one bgen chunk for this chromosome"),
  make_option("--sample",      type = "character", help = "its .sample file"),
  make_option("--variants",    type = "character", help = "variant list to include: tsv with chr,pos,a0,a1 (e.g. HapMap3+)"),
  make_option("--out",         type = "character", help = "output per-chr bundle .rds"),
  make_option("--n_ld",        type = "integer",   default = 10000L, help = "individuals subsampled for LD [default %default]"),
  make_option("--ld_ids",      type = "character", default = NULL, help = "optional file of sample ids to use (overrides --n_ld)"),
  make_option("--genetic_map", type = "character", default = NULL, help = "optional genetic-map dir for snp_asGeneticPos (cM window)"),
  make_option("--window_cM",   type = "double",    default = 3,    help = "LD window in cM when a genetic map is given [default %default]"),
  make_option("--window_kb",   type = "double",    default = 3000, help = "LD window in kb when no genetic map [default %default]"),
  make_option("--ncores",      type = "integer",   default = 1L)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("bgen","sample","variants","out")) if (is.null(opt[[r]])) stop("Required: --bgen --sample --variants --out")
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

chr_to_int <- function(x) {
  x <- sub("^CHR", "", toupper(as.character(x)))
  x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x == "XY"] <- "25"; x[x %in% c("MT","M")] <- "26"
  suppressWarnings(as.integer(x))
}
allele_key <- function(a, b) { a <- toupper(a); b <- toupper(b); ifelse(a < b, paste(a,b,sep="/"), paste(b,a,sep="/")) }

## ---- variant list to include -------------------------------------------
vars <- fread(opt$variants)
setnames(vars, tolower(names(vars)))
if (!all(c("chr","pos","a0","a1") %in% names(vars))) stop("--variants needs columns chr,pos,a0,a1")
vars_key <- paste(chr_to_int(vars$chr), vars$pos, allele_key(vars$a1, vars$a0), sep = ":")

## ---- .bgi: read + rewrite chromosome to bare integer (see score_chunk.R) --
bgi <- paste0(opt$bgen, ".bgi")
con <- dbConnect(RSQLite::SQLite(), bgi)
vt <- as.data.table(dbGetQuery(con, "SELECT chromosome, position, allele1, allele2 FROM Variant"))
for (ch in unique(vt$chromosome))
  dbExecute(con, "UPDATE Variant SET chromosome = :n WHERE chromosome = :o",
            params = list(n = as.character(chr_to_int(ch)), o = ch))
dbDisconnect(con)
bgi_key <- paste(chr_to_int(vt$chromosome), vt$position, allele_key(vt$allele1, vt$allele2), sep = ":")
bgi_ids <- paste(chr_to_int(vt$chromosome), vt$position, vt$allele1, vt$allele2, sep = "_")
want <- unique(bgi_ids[bgi_key %in% vars_key])
log_msg("chr variants in list present in chunk: %d", length(want))
if (length(want) == 0) stop("no listed variants present in this chunk")

## ---- individual subsample ----------------------------------------------
samp <- fread(opt$sample, header = FALSE, skip = 2)
ids  <- as.character(if (ncol(samp) >= 2) samp[[2]] else samp[[1]])
n_tot <- length(ids)
if (!is.null(opt$ld_ids)) {
  use <- readLines(opt$ld_ids)
  ind_row <- which(ids %in% use)
} else {
  set.seed(1)
  ind_row <- sort(sample(n_tot, min(opt$n_ld, n_tot)))
}
log_msg("using %d / %d individuals for LD", length(ind_row), n_tot)

## ---- read dosages + compute correlation ---------------------------------
backing <- tempfile(tmpdir = getwd())
rds <- snp_readBGEN(bgenfiles = opt$bgen, backingfile = backing,
                    list_snp_id = list(want), ind_row = ind_row, ncores = opt$ncores)
obj <- snp_attach(rds)
G <- obj$genotypes
map <- obj$map
## Drop variants monomorphic in the LD subsample: snp_cor returns NaN for
## zero-variance columns, which would poison the SFBM / snp_ldpred2_auto.
cs   <- bigstatsr::big_colstats(G, ncores = opt$ncores)
keep <- which(cs$var > 0)
if (length(keep) < ncol(G)) log_msg("dropping %d monomorphic-in-subsample variant(s)", ncol(G) - length(keep))
if (length(keep) == 0) stop("all variants monomorphic in the LD subsample")
POS <- map$physical.pos[keep]
## bigsnpr::snp_cor `size` is a window of "1000 * size" positions in the units of
## infos.pos. So with bp positions, size = window in kb; with cM positions,
## size = window_cM / 1000 (cf. the LDpred2 vignette's `size = 3 / 1000`).
if (!is.null(opt$genetic_map)) {
  POS2 <- snp_asGeneticPos(chr_to_int(map$chromosome[keep]), POS, dir = opt$genetic_map, ncores = opt$ncores)
  corr <- snp_cor(G, ind.col = keep, infos.pos = POS2, size = opt$window_cM / 1000, ncores = opt$ncores)
} else {
  corr <- snp_cor(G, ind.col = keep, infos.pos = POS, size = opt$window_kb, ncores = opt$ncores)  # POS bp, size kb
}
ld <- Matrix::colSums(corr^2)

out_map <- data.frame(chr = chr_to_int(map$chromosome[keep]), pos = map$physical.pos[keep],
                      a0 = map$allele1[keep], a1 = map$allele2[keep], stringsAsFactors = FALSE)
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
saveRDS(list(corr = corr, map = out_map, ld = ld), opt$out)
log_msg("wrote per-chr LD bundle: %d variants -> %s", nrow(out_map), opt$out)
