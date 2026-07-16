#!/usr/bin/env Rscript

## ldpred2_infer.R
##
## LDpred2-auto weight INFERENCE. Given GWAS summary statistics and an LD
## reference (per-chromosome correlation "parts"), infers per-SNP posterior
## effect sizes with bigsnpr::snp_ldpred2_auto and writes them as a harmonised
## weight file that scripts/score_chunk.R can score:
##     bgen_id  chr  pos  a1_effect  a0  weight
##
## Inference half of the LDpred2 pipeline (scoring half = score_chunk.R +
## gather_scores.R). Only summary statistics and the LD reference are used here.
##
## LD reference = one or more PART files (from build_ld_reference.R, or a
## precomputed reference converted to this format). Each part is an .rds holding
##   list(corr = <sparse symmetric correlation, dsCMatrix>,
##        map  = data.frame(chr, pos, a0, a1)   # rows match corr,
##        ld   = <numeric LD scores, one per variant>)
## typically one part per chromosome. Correlation is block-diagonal across parts.
##
## Portability note: parts store PLAIN sparse matrices (portable across machines).
## The on-disk SFBM that snp_ldpred2_auto needs is (re)built HERE, in this task's
## working dir, after matching + subsetting to the sumstats -- so no absolute
## backingfile path ever has to survive being moved between WDL tasks.
##
## Canonical workflow: per-part snp_match -> subset+assemble SFBM aligned to
## df_beta -> snp_ldsc (h2) -> snp_ldpred2_auto (multi-chain) -> filter diverged
## chains -> average beta_est -> weights.

suppressMessages({
  library(optparse)
  library(data.table)
  library(bigsnpr)
  library(bigsparser)
  library(Matrix)
})
log_msg <- function(...) cat(sprintf("[ldpred2_infer] %s\n", sprintf(...)), file = stderr())

option_list <- list(
  make_option("--sumstats",   type = "character", help = "GWAS summary statistics (tsv/csv, optionally .gz)"),
  make_option("--ld_parts",   type = "character", help = "per-chr LD part .rds files: manifest (one path/line) or comma list; or a single bundle .rds with $parts"),
  make_option("--out",        type = "character", help = "output harmonised weight file (.tsv[.gz])"),
  make_option("--trait_name", type = "character", default = "trait"),
  make_option("--chr_col",    type = "character", default = "chr"),
  make_option("--pos_col",    type = "character", default = "pos"),
  make_option("--a1_col",     type = "character", default = "a1", help = "effect allele column"),
  make_option("--a0_col",     type = "character", default = "a0", help = "other allele column"),
  make_option("--beta_col",   type = "character", default = "beta"),
  make_option("--se_col",     type = "character", default = "beta_se"),
  make_option("--n_col",      type = "character", default = "n_eff", help = "per-variant effective N column"),
  make_option("--freq_col",   type = "character", default = NULL, help = "a1 allele-freq column (enables sd-based QC)"),
  make_option("--info_col",   type = "character", default = NULL, help = "imputation INFO column"),
  make_option("--n_eff",      type = "double",    default = NA,   help = "constant effective N (if --n_col absent)"),
  make_option("--info_thr",   type = "double",    default = 0.7),
  make_option("--sd_qc",      type = "double",    default = 0.10),
  make_option("--h2_init",    type = "double",    default = NA, help = "override the LDSC h2 seed (use when LDSC is unreliable, e.g. few variants)"),
  make_option("--n_chains",   type = "integer",   default = 30L),
  make_option("--burn_in",    type = "integer",   default = 500L),
  make_option("--num_iter",   type = "integer",   default = 500L),
  make_option("--ncores",     type = "integer",   default = 1L),
  make_option("--report",     type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))
for (r in c("sumstats","ld_parts","out")) if (is.null(opt[[r]])) stop("Required: --sumstats --ld_parts --out")
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

chr_to_int <- function(x) {
  x <- sub("^CHR", "", toupper(as.character(x)))
  x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x == "XY"] <- "25"; x[x %in% c("MT","M")] <- "26"
  suppressWarnings(as.integer(x))
}

## ---- 1. read + shape summary statistics ---------------------------------
log_msg("reading sumstats: %s", opt$sumstats)
ss <- fread(opt$sumstats)
need <- c(opt$chr_col, opt$pos_col, opt$a1_col, opt$a0_col, opt$beta_col, opt$se_col)
miss <- setdiff(need, names(ss)); if (length(miss)) stop(sprintf("sumstats missing columns: %s", paste(miss, collapse=", ")))
sumstats <- data.frame(
  chr = chr_to_int(ss[[opt$chr_col]]), pos = as.integer(ss[[opt$pos_col]]),
  a1 = toupper(as.character(ss[[opt$a1_col]])), a0 = toupper(as.character(ss[[opt$a0_col]])),
  beta = as.numeric(ss[[opt$beta_col]]), beta_se = as.numeric(ss[[opt$se_col]]),
  stringsAsFactors = FALSE)
if (!is.null(opt$n_col) && opt$n_col %in% names(ss)) sumstats$n_eff <- as.numeric(ss[[opt$n_col]]) else
  if (is.finite(opt$n_eff)) sumstats$n_eff <- opt$n_eff else stop("provide --n_col or --n_eff")
if (!is.null(opt$freq_col) && opt$freq_col %in% names(ss)) sumstats$freq <- as.numeric(ss[[opt$freq_col]])
if (!is.null(opt$info_col) && opt$info_col %in% names(ss)) sumstats$info <- as.numeric(ss[[opt$info_col]])
n0 <- nrow(sumstats)
sumstats <- sumstats[is.finite(sumstats$chr) & is.finite(sumstats$pos) &
                       is.finite(sumstats$beta) & is.finite(sumstats$beta_se) & sumstats$beta_se > 0, , drop = FALSE]
if (!is.null(sumstats$info)) sumstats <- sumstats[is.na(sumstats$info) | sumstats$info >= opt$info_thr, , drop = FALSE]
log_msg("sumstats: %d rows (%d after basic filters)", n0, nrow(sumstats))

## ---- 2. resolve LD parts ------------------------------------------------
resolve_parts <- function(spec) {
  if (file.exists(spec) && grepl("\\.rds$", spec, ignore.case = TRUE) && !grepl(",", spec)) {
    b <- readRDS(spec)
    if (!is.null(b$parts)) return(b$parts)            # a bundle of in-memory parts
    if (!is.null(b$corr) && !is.null(b$map)) return(list(b))  # a single part
    stop("unrecognised LD .rds (need $parts, or $corr+$map)")
  }
  paths <- if (file.exists(spec) && !grepl(",", spec)) readLines(spec) else strsplit(spec, ",")[[1]]
  paths <- trimws(paths); paths <- paths[nchar(paths) > 0 & file.exists(paths)]
  if (length(paths) == 0) stop("no LD part files found")
  lapply(paths, readRDS)
}
parts <- resolve_parts(opt$ld_parts)
log_msg("LD reference: %d part(s)", length(parts))

## ---- 3. per-part: match -> (canonical sd-QC) -> subset -> assemble SFBM ---
## QC is applied WITHIN each block, before the block's correlation is appended,
## so the SFBM is built already-clean and aligned to df_beta (we never subset an
## SFBM after the fact). Canonical bigsnpr sd-QC: the GWAS-implied allele SD
## sd_ss = 2/sqrt(n_eff*se^2 + beta^2) must agree with sd_af = sqrt(2*f*(1-f)).
tmp_sfbm <- tempfile(tmpdir = getwd())
corr <- NULL
df_list <- list(); info_list <- list(); ld_keep <- numeric(0)
n_matched <- 0L; n_qc_drop <- 0L; did_qc <- FALSE; n_ref_total <- 0L
for (k in seq_along(parts)) {
  cc <- parts[[k]]$corr
  m  <- as.data.frame(parts[[k]]$map); ldv <- as.numeric(parts[[k]]$ld)
  m$chr <- chr_to_int(m$chr)
  n_ref_total <- n_ref_total + nrow(m)          # full LD-reference variant count M (for LDSC)
  mm <- tryCatch(snp_match(sumstats, m, join_by_pos = TRUE, match.min.prop = 0),
                 error = function(e) NULL)              # strand_flip default TRUE (external sumstats)
  if (is.null(mm) || nrow(mm) == 0) next
  n_matched <- n_matched + nrow(mm)
  ok <- rep(TRUE, nrow(mm))
  if (!is.null(mm$freq)) {
    did_qc <- TRUE
    sd_af <- sqrt(2 * mm$freq * (1 - mm$freq))
    sd_ss <- 2 / sqrt(mm$n_eff * mm$beta_se^2 + mm$beta^2)
    ok <- !(sd_ss < 0.5 * sd_af | sd_ss > sd_af + opt$sd_qc | sd_ss < 0.1 | sd_af < 0.05)
    ok[is.na(ok)] <- FALSE
    n_qc_drop <- n_qc_drop + sum(!ok)
  }
  if (!any(ok)) next
  mm <- mm[ok, , drop = FALSE]
  loc <- mm[["_NUM_ID_"]]                                # indices into m / cc, in mm order
  cc_sub <- cc[loc, loc, drop = FALSE]
  cc_sub <- Matrix::forceSymmetric(cc_sub)               # robust symmetric coercion -> dsCMatrix
  if (is.null(corr)) corr <- as_SFBM(cc_sub, tmp_sfbm, compact = TRUE)
  else               corr$add_columns(cc_sub, nrow(corr))
  df_list[[length(df_list) + 1L]]   <- data.frame(beta = mm$beta, beta_se = mm$beta_se, n_eff = mm$n_eff)
  info_list[[length(info_list) + 1L]] <- data.frame(chr = mm$chr, pos = mm$pos, a1 = mm$a1, a0 = mm$a0)
  ld_keep <- c(ld_keep, ldv[loc])
}
if (is.null(corr)) stop("no sumstats variants matched the LD reference")
df_beta  <- do.call(rbind, df_list)
info_snp <- do.call(rbind, info_list)
if (did_qc) { log_msg("sd-based QC: dropped %d / %d matched variants", n_qc_drop, n_matched)
} else { log_msg("sd-based QC skipped (no --freq_col)") }
log_msg("assembled %d variants across %d block(s) (SFBM ncol=%d)", nrow(df_beta), length(df_list), ncol(corr))
if (nrow(df_beta) != ncol(corr)) stop("internal: df_beta rows != corr columns")

## ---- 5. heritability (LDSC) ---------------------------------------------
chi2 <- (df_beta$beta / df_beta$beta_se)^2
ldsc <- snp_ldsc(ld_score = ld_keep, ld_size = n_ref_total, chi2 = chi2,
                 sample_size = df_beta$n_eff, ncores = opt$ncores)
h2_est <- ldsc[["h2"]]
log_msg("LDSC h2 = %.4f (intercept %.3f)", h2_est, ldsc[["int"]])
if (!is.finite(h2_est) || h2_est <= 0) { log_msg("h2 non-positive/NA -> fallback 0.1"); h2_est <- 0.1 }
if (is.finite(opt$h2_init) && opt$h2_init > 0) {
  log_msg("using provided --h2_init %.4f (overriding LDSC %.4f)", opt$h2_init, h2_est); h2_est <- opt$h2_init
}

## ---- 6. LDpred2-auto ----------------------------------------------------
set.seed(1)
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = opt$n_chains),
  burn_in = opt$burn_in, num_iter = opt$num_iter,
  ncores = opt$ncores, allow_jump_sign = FALSE, shrink_corr = 0.95)
rng  <- sapply(multi_auto, function(a) diff(range(a$corr_est, na.rm = TRUE)))
keep <- which(rng > (0.95 * stats::quantile(rng, 0.95, na.rm = TRUE)))
if (length(keep) == 0) { log_msg("WARNING: no chains passed filter; keeping all"); keep <- seq_along(multi_auto) }
log_msg("LDpred2-auto: kept %d / %d chains", length(keep), length(multi_auto))
## do.call(cbind, ...) keeps it a matrix even when a single chain is kept
## (sapply would drop to a vector and break rowMeans).
beta_auto <- rowMeans(do.call(cbind, lapply(multi_auto[keep], function(a) a$beta_est)))
p_est    <- mean(sapply(multi_auto[keep], function(a) a$p_est),  na.rm = TRUE)
h2_final <- mean(sapply(multi_auto[keep], function(a) a$h2_est), na.rm = TRUE)

## ---- 7. write harmonised weights ----------------------------------------
kv <- is.finite(beta_auto) & beta_auto != 0
out <- data.table(
  bgen_id   = paste(info_snp$chr, info_snp$pos, info_snp$a0, info_snp$a1, sep = "_")[kv],
  chr = info_snp$chr[kv], pos = info_snp$pos[kv],
  a1_effect = info_snp$a1[kv], a0 = info_snp$a0[kv], weight = beta_auto[kv])
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
if (grepl("\\.gz$", opt$out)) {
  tmp <- sub("\\.gz$", "", opt$out); fwrite(out, tmp, sep = "\t", quote = FALSE)
  system2("gzip", c("-f", shQuote(tmp)))     # image data.table lacks zlib
} else fwrite(out, opt$out, sep = "\t", quote = FALSE)
log_msg("wrote %d weights -> %s", nrow(out), opt$out)

if (!is.null(opt$report)) {
  fwrite(data.table(trait = opt$trait_name, n_sumstats = n0, n_matched = nrow(df_beta),
                    n_weights = nrow(out), h2_ldsc = h2_est, h2_auto = h2_final,
                    p_est = p_est, n_chains_kept = length(keep)),
         opt$report, sep = "\t", quote = FALSE)
}
log_msg("done: %s (h2=%.4f, p=%.3g, %d weights)", opt$trait_name, h2_final, p_est, nrow(out))
