# AANAT7223: GSE23179 leave-one-array-out sensitivity analysis
#
# Purpose:
#   Test whether the 15 primary GSE23179 dose-vs-time-matched-control contrasts
#   are materially influenced by any single array, including the 13 arrays
#   flagged by the exploratory RLE/NUSE screen.
#
# Important design choice:
#   The common 72-array raw-CEL/RMA expression matrix is held fixed, and limma
#   is refitted 72 times after omitting one array at a time. This isolates
#   sample-level influence on the limma estimates while avoiding a different
#   quantile-normalization reference in every iteration. It is NOT described as
#   an RMA re-normalization sensitivity analysis.
#
# Inputs:
#   - AANAT7223_AllFeature_EffectSize_95CI_Audit_v1.zip
#   - Optional GSE23179_RLE_NUSE_sample_level_review.csv
#
# The script downloads GSE23179 raw CEL files and Series Matrix metadata from GEO
# if they are not already cached. All 72 arrays remain in the primary analysis.

options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20261109)

# ------------------------- User-configurable paths -------------------------
INPUT_ZIP <- "AANAT7223_AllFeature_EffectSize_95CI_Audit_v1.zip"
QC_CSV <- "GSE23179_RLE_NUSE_sample_level_review.csv"
OUT <- "AANAT7223_GSE23179_LOAO_v2"
RAW_DIR <- file.path(OUT, "rawCEL")
MATRIX_DIR <- file.path(OUT, "series_matrix")
EXPR_DIR <- file.path(OUT, "RMA_expression")
RESULT_DIR <- file.path(OUT, "results")
for (d in c(OUT, RAW_DIR, MATRIX_DIR, EXPR_DIR, RESULT_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ Packages -----------------------------------
# Use the packages already present in the project runtime. The cached RMA matrix
# and Series Matrix are sufficient for LOAO; GEOquery/affy are not required.
for (p in c("data.table", "limma")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Required package is not installed in this runtime: ", p,
         ". Install it in the project's supported Bioconductor environment.")
}
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

# ------------------------------ GEO inputs ---------------------------------
# Read the sample manifest directly from the cached GEO Series Matrix. This
# avoids a GEOquery dependency and does not alter any sample metadata.
SERIES_MATRIX <- file.path(
  "AANAT7223_Ranked_GSEA_Camera_v2", "rawCEL", "GSE23179_matrix",
  "GSE23179_series_matrix.txt.gz"
)
if (!file.exists(SERIES_MATRIX))
  stop("Missing cached GSE23179 Series Matrix: ", SERIES_MATRIX)
con <- gzfile(SERIES_MATRIX, open = "rt")
series_lines <- readLines(con, warn = FALSE)
close(con)

extract_geo_sample_field <- function(key) {
  prefix <- paste0(key, "\t")
  hit <- series_lines[startsWith(series_lines, prefix)]
  if (length(hit) != 1L)
    stop("Expected exactly one ", key, " field in Series Matrix; found ", length(hit))
  values <- strsplit(hit, "\t", fixed = TRUE)[[1]]
  values <- values[-1L]
  values <- sub('^"', "", values)
  values <- sub('"$', "", values)
  values
}
manifest <- data.table(
  GSM = toupper(extract_geo_sample_field("!Sample_geo_accession")),
  title = extract_geo_sample_field("!Sample_title")
)
if (nrow(manifest) != 72L || uniqueN(manifest$GSM) != 72L ||
    anyNA(manifest$title) || any(!nzchar(manifest$title)))
  stop("Series Matrix metadata must contain exactly 72 unique GSM accessions and valid titles.")

# Reuse the RMA matrix produced and audited in the four-dataset workflow.
expr_file <- file.path(
  "AANAT7223_Ranked_GSEA_Camera_v2", "RMA_expression",
  "GSE23179_RMA_expression.rds"
)
if (!file.exists(expr_file))
  stop("Missing cached RMA expression matrix: ", expr_file,
       ". Do not silently substitute a different normalization workflow.")
expr <- readRDS(expr_file)
if (nrow(expr) != 31099L || ncol(expr) != 72L)
  stop("Expected 31,099 probe sets x 72 arrays; observed ",
       nrow(expr), " x ", ncol(expr), ".")
if (is.null(colnames(expr)) || !setequal(colnames(expr), manifest$GSM))
  stop("RMA matrix column names do not match the 72 GSM accessions from Series Matrix.")
expr <- expr[, manifest$GSM, drop = FALSE]

# -------------------------- Exposure group design --------------------------
# The GEO title uses EB_07 for 0.7 ppm and EB_0 for control.
dose_raw <- sub(".*EB_([^_]+)_.*", "\\1", manifest$title)
dose_code <- c("0" = "0", "07" = "0p7", "0.7" = "0p7",
               "2" = "2", "6" = "6", "10" = "10", "15" = "15")[dose_raw]
duration_code <- ifelse(grepl("_5d_", manifest$title), "Td5",
                 ifelse(grepl("_4w_", manifest$title), "Tw4",
                 ifelse(grepl("_13w_", manifest$title), "Tw13", NA_character_)))
if (anyNA(dose_code) || anyNA(duration_code))
  stop("Could not parse dose/duration from all GEO sample titles.")
group <- paste0(duration_code, "_D", dose_code)
group_levels <- unlist(lapply(c("Td5", "Tw4", "Tw13"), function(d)
  paste0(d, "_D", c("0", "0p7", "2", "6", "10", "15"))),
  use.names = FALSE)
group <- factor(group, levels = group_levels)
group_counts <- as.data.table(table(group))
setnames(group_counts, c("group", "n"))
if (nrow(group_counts) != 18L || any(group_counts$n != 4L))
  stop("Expected 18 dose-duration groups with 4 arrays each.")

# Construct contrasts with the exact labels used in the V23 all-feature archive.
dose_contrasts <- c("0p7", "2", "6", "10", "15")
contrast_defs <- rbindlist(lapply(c("Td5", "Tw4", "Tw13"), function(dp) {
  data.table(
    Contrast = paste0(dp, "_D", dose_contrasts, ".", dp, "_D0"),
    Exposed_group = paste0(dp, "_D", dose_contrasts),
    Control_group = paste0(dp, "_D0")
  )
}))
cm <- matrix(0, nrow = length(group_levels), ncol = nrow(contrast_defs),
             dimnames = list(group_levels, contrast_defs$Contrast))
for (j in seq_len(nrow(contrast_defs))) {
  cm[contrast_defs$Exposed_group[j], j] <- 1
  cm[contrast_defs$Control_group[j], j] <- -1
}

# --------------------------- Primary model fit ------------------------------
design_primary <- model.matrix(~ 0 + group)
colnames(design_primary) <- group_levels
fit_primary <- limma::lmFit(expr, design_primary)
fit_primary <- limma::eBayes(limma::contrasts.fit(fit_primary, cm),
                             trend = TRUE, robust = TRUE)
base_tabs <- vector("list", nrow(contrast_defs))
names(base_tabs) <- contrast_defs$Contrast
for (j in seq_len(nrow(contrast_defs))) {
  tt <- limma::topTable(fit_primary, coef = j, number = Inf,
                        adjust.method = "BH", sort.by = "none")
  base_tabs[[j]] <- data.table(
    Probe_set = rownames(tt),
    log2FC = tt$logFC,
    Moderated_t = tt$t,
    BH_FDR = tt$adj.P.Val
  )
}

# ---------------------- Validate against locked V23 ------------------------
stat_dir <- file.path(OUT, "all_feature_statistics")
dir.create(stat_dir, showWarnings = FALSE, recursive = TRUE)
unzip(INPUT_ZIP, exdir = stat_dir)
stat_file <- list.files(stat_dir,
  pattern = "GSE23179_AllFeature_Limma_Conservative95CI\\.csv\\.gz$",
  recursive = TRUE, full.names = TRUE)
if (length(stat_file) != 1L) stop("Could not uniquely locate the GSE23179 V23 all-feature table.")
locked <- data.table::fread(stat_file)
locked[, Probe_set := as.character(Probe_set)]
validation <- rbindlist(lapply(names(base_tabs), function(cn) {
  b <- base_tabs[[cn]]
  e <- locked[Contrast == cn, .(Probe_set, log2FC, Moderated_t, BH_FDR)]
  if (nrow(e) != 31099L) stop("Unexpected locked V23 row count for ", cn)
  idx <- match(b$Probe_set, e$Probe_set)
  if (anyNA(idx) || length(idx) != 31099L) stop("Feature ID mismatch for ", cn)
  e_sub <- e[idx]
  data.table(
    Contrast = cn,
    N_features = nrow(b),
    Spearman_t_rerun_vs_V23 = cor(b$Moderated_t, e_sub$Moderated_t,
                                  method = "spearman", use = "complete.obs"),
    Max_abs_log2FC_difference = max(abs(b$log2FC - e_sub$log2FC), na.rm = TRUE),
    FDR05_rerun = sum(b$BH_FDR < 0.05, na.rm = TRUE),
    FDR05_V23 = sum(e_sub$BH_FDR < 0.05, na.rm = TRUE)
  )
}))
fwrite(validation, file.path(RESULT_DIR, "GSE23179_LOAO_baseline_validation.csv"))
if (any(validation$Spearman_t_rerun_vs_V23 < 0.99999) ||
    any(validation$Max_abs_log2FC_difference > 1e-5) ||
    any(validation$FDR05_rerun != validation$FDR05_V23)) {
  stop("The reconstructed baseline does not sufficiently reproduce the locked V23 results. " ,
       "Review GSE23179_LOAO_baseline_validation.csv before interpreting LOAO.")
}

# ---------------------------- QC flag lookup -------------------------------
qc <- NULL
if (file.exists(QC_CSV)) {
  qc <- fread(QC_CSV)
  if (!"GSM" %in% names(qc) || !"robust_screen_flag" %in% names(qc))
    stop("QC_CSV must contain GSM and robust_screen_flag columns.")
  qc[, GSM := toupper(GSM)]
} else {
  qc <- data.table(GSM = manifest$GSM, robust_screen_flag = NA)
}
qc[, robust_screen_flag := as.logical(robust_screen_flag)]
qc_flag <- qc$robust_screen_flag[match(manifest$GSM, qc$GSM)]

# ----------------------- Leave-one-array-out fits --------------------------
loao_rows <- vector("list", ncol(expr) * nrow(contrast_defs))
out_i <- 0L
for (i in seq_len(ncol(expr))) {
  message("LOAO: omitting ", i, "/", ncol(expr), " ", manifest$GSM[i])
  keep <- setdiff(seq_len(ncol(expr)), i)
  group_i <- factor(as.character(group[keep]), levels = group_levels)
  design_i <- model.matrix(~ 0 + group_i)
  colnames(design_i) <- group_levels
  fit_i <- limma::lmFit(expr[, keep, drop = FALSE], design_i)
  fit_i <- limma::eBayes(limma::contrasts.fit(fit_i, cm),
                         trend = TRUE, robust = TRUE)

  for (j in seq_len(nrow(contrast_defs))) {
    cn <- contrast_defs$Contrast[j]
    tt <- limma::topTable(fit_i, coef = j, number = Inf,
                          adjust.method = "BH", sort.by = "none")
    loo <- data.table(
      Probe_set = rownames(tt),
      log2FC = tt$logFC,
      Moderated_t = tt$t,
      BH_FDR = tt$adj.P.Val
    )
    b <- base_tabs[[cn]]
    if (!identical(loo$Probe_set, b$Probe_set))
      loo <- loo[match(b$Probe_set, loo$Probe_set)]
    if (anyNA(loo$Probe_set)) stop("Feature matching failed in LOAO for ", cn)

    delta <- abs(loo$log2FC - b$log2FC)
    base_sig <- which(b$BH_FDR < 0.05)
    loo_sig <- which(loo$BH_FDR < 0.05)
    union_sig <- union(base_sig, loo_sig)
    jaccard <- if (!length(union_sig)) NA_real_ else length(intersect(base_sig, loo_sig)) / length(union_sig)
    top500 <- order(abs(b$Moderated_t), decreasing = TRUE)[seq_len(500L)]
    sign_switch <- if (!length(base_sig)) 0L else
      sum(sign(b$log2FC[base_sig]) != sign(loo$log2FC[base_sig]), na.rm = TRUE)

    out_i <- out_i + 1L
    loao_rows[[out_i]] <- data.table(
      Omitted_GSM = manifest$GSM[i],
      Omitted_group = as.character(group[i]),
      Omitted_PLM_flag = ifelse(is.na(qc_flag[i]), NA, qc_flag[i]),
      Contrast = cn,
      Omitted_sample_role = if (as.character(group[i]) == contrast_defs$Exposed_group[j]) "exposed" else
                            if (as.character(group[i]) == contrast_defs$Control_group[j]) "control" else "other",
      Spearman_t_vs_primary = cor(b$Moderated_t, loo$Moderated_t,
                                  method = "spearman", use = "complete.obs"),
      Pearson_log2FC_vs_primary = cor(b$log2FC, loo$log2FC,
                                      method = "pearson", use = "complete.obs"),
      Median_abs_delta_log2FC = median(delta, na.rm = TRUE),
      P95_abs_delta_log2FC = as.numeric(quantile(delta, 0.95, na.rm = TRUE)),
      Max_abs_delta_log2FC = max(delta, na.rm = TRUE),
      Primary_FDR05_features = length(base_sig),
      LOAO_FDR05_features = length(loo_sig),
      Delta_FDR05_features = length(loo_sig) - length(base_sig),
      FDR05_Jaccard = jaccard,
      Primary_FDR05_direction_switches = sign_switch,
      Top500_primary_direction_concordance = mean(
        sign(b$log2FC[top500]) == sign(loo$log2FC[top500]), na.rm = TRUE
      )
    )
  }
}
loao <- rbindlist(loao_rows, fill = TRUE)
fwrite(loao, file.path(RESULT_DIR, "GSE23179_LOAO_all_72_arrays_x_15_contrasts.csv"))

# ------------------------- Contrast-level summary --------------------------
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
summary <- loao[, .(
  N_LOAO_fits = .N,
  Min_Spearman_t = min(Spearman_t_vs_primary, na.rm = TRUE),
  Median_Spearman_t = median(Spearman_t_vs_primary, na.rm = TRUE),
  Min_Pearson_log2FC = min(Pearson_log2FC_vs_primary, na.rm = TRUE),
  Max_P95_abs_delta_log2FC = max(P95_abs_delta_log2FC, na.rm = TRUE),
  Max_abs_delta_log2FC = max(Max_abs_delta_log2FC, na.rm = TRUE),
  Min_FDR05_Jaccard = safe_min(FDR05_Jaccard),
  Min_LOAO_FDR05_features = min(LOAO_FDR05_features, na.rm = TRUE),
  Max_LOAO_FDR05_features = max(LOAO_FDR05_features, na.rm = TRUE),
  Max_abs_FDR05_count_change = safe_max(abs(Delta_FDR05_features)),
  Max_primary_sign_switches = max(Primary_FDR05_direction_switches, na.rm = TRUE),
  Min_Top500_direction_concordance = min(Top500_primary_direction_concordance, na.rm = TRUE)
), by = Contrast]
setorder(summary, Contrast)
fwrite(summary, file.path(RESULT_DIR, "GSE23179_LOAO_by_contrast_summary.csv"))

flag_summary <- loao[, .(
  N_omissions = .N,
  Median_Spearman_t = median(Spearman_t_vs_primary, na.rm = TRUE),
  Min_Spearman_t = min(Spearman_t_vs_primary, na.rm = TRUE),
  Median_FDR05_Jaccard = safe_median(FDR05_Jaccard),
  Min_FDR05_Jaccard = safe_min(FDR05_Jaccard),
  Max_abs_FDR05_count_change = max(abs(Delta_FDR05_features), na.rm = TRUE)
), by = .(Omitted_PLM_flag)]
fwrite(flag_summary, file.path(RESULT_DIR, "GSE23179_LOAO_by_PLM_flag_status.csv"))

# ------------------------------- Provenance --------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
writeLines(c(
  "GSE23179 leave-one-array-out limma sensitivity analysis.",
  "All 72 arrays were retained in the primary model.",
  "Sensitivity: 72 separate refits, each omitting one array from the common 72-array raw-CEL/RMA expression matrix.",
  "RMA normalization was held fixed across LOAO refits; this isolates sample-level influence on limma estimates.",
  "All 15 time-matched dose-vs-control contrasts were refitted in each omission.",
  "Baseline reproduction against the locked V23 all-feature table was required before LOAO interpretation.",
  "FDR was recalculated across the full 31,099-feature universe for every contrast and every LOAO fit.",
  "The 13 RLE/NUSE flags are exploratory review indicators, not automatic exclusion criteria."
), file.path(OUT, "README_analysis_provenance.txt"))

if (nzchar(Sys.which("zip"))) {
  old <- getwd(); setwd(OUT)
  zip("../AANAT7223_GSE23179_LOAO_v2_results.zip",
      files = c("results", "sessionInfo.txt", "README_analysis_provenance.txt"),
      flags = "-r")
  setwd(old)
}
message("Completed GSE23179 leave-one-array-out sensitivity analysis.")
message("Review baseline validation before interpreting the LOAO tables.")