# GSE21477 source-ID-blocked reanalysis, version 2
# Annals of Anatomy major revision
#
# Purpose
#   1) Rebuild the 39-sample manifest from GEO metadata (never infer groups by column position).
#   2) Reprocess raw CEL files with Affymetrix RMA.
#   3) Treat ES1, ES2 and ES3 as distinct experimental blocks.
#   4) Primary inferential analyses: source-ID fixed-effect limma models fitted separately
#      within ES1, ES2 and ES3, using the common 39-array RMA matrix.
#   5) Sensitivity analyses: (a) all-block source-ID-nested model on the common RMA matrix;
#      (b) historical all-block unblocked model; (c) within-block RMA followed by
#      block-specific source-ID models.
#   6) Generate chip-level QC, including RLE/NUSE when affyPLM is available.
#
# Important interpretation rules
#   - HNEC source identifiers are used only within their GEO experimental block.
#   - Do not assume source IDs are globally interchangeable across ES1/ES2/ES3.
#   - ES1, ES2 and ES3 have different treatment/control structures; their contrasts
#     must remain block-specific.
#   - Micromolar in-vitro concentrations must not be converted to inhaled ppm.
#   - This script must be executed in a complete R/Bioconductor environment. No results
#     are claimed until the script has been run and its outputs independently checked.

options(stringsAsFactors = FALSE, timeout = 600)

# ---------- 0. Packages ----------
cran_pkgs <- c("data.table", "matrixStats", "R.utils")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
bioc_pkgs <- c("GEOquery", "Biobase", "affy", "limma", "hgu133plus2cdf")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
if (!requireNamespace("affyPLM", quietly = TRUE)) {
  try(BiocManager::install("affyPLM", ask = FALSE, update = FALSE), silent = TRUE)
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(affy)
  library(limma)
  library(data.table)
  library(matrixStats)
  library(R.utils)
})

OUT <- "GSE21477_sourceID_blocked_limma_v2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
for (d in c("rawCEL", "global_RMA", "within_block_RMA", "qc", "metadata", "logs")) {
  dir.create(file.path(OUT, d), showWarnings = FALSE, recursive = TRUE)
}

# ---------- 1. Metadata from GEO ----------
gse_list <- getGEO("GSE21477", GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gse_list) != 1L) stop("Expected one GPL570 SeriesMatrix object.")
gse <- gse_list[[1]]
pd <- pData(gse)
if (!all(c("title", "source_name_ch1") %in% names(pd))) {
  stop("GEO metadata is missing title or source_name_ch1.")
}
md <- data.frame(
  GSM = rownames(pd),
  title = as.character(pd$title),
  source = as.character(pd$source_name_ch1),
  stringsAsFactors = FALSE
)
md$block <- sub("^(ES[123]).*", "\\1", md$source)
md$source_ID <- sub(".*HNEC_#([0-9]+).*", "\\1", md$source)
md$bio_rep <- sub(".*biological rep([0-9]+).*", "\\1", md$title, ignore.case = TRUE)

classify_condition <- function(title, block) {
  tt <- tolower(title)
  if (block == "ES1") {
    if (grepl("50.*fa for 2 hours", tt)) return("50uM_2h")
    if (grepl("50.*fa for 4 hours", tt)) return("50uM_4h")
    if (grepl("100.*fa for 2 hours", tt)) return("100uM_2h")
    if (grepl("control", tt)) return("control_shared")
  }
  if (block == "ES2") {
    if (grepl("100.*fa for 4 hours", tt)) return("100uM_4h")
    if (grepl("100.*fa for 24 hours", tt)) return("100uM_24h")
    if (grepl("200.*fa for 4 hours", tt)) return("200uM_4h")
    if (grepl("200.*fa for 24 hours", tt)) return("200uM_24h")
    if (grepl("control group after 4 hours", tt)) return("control_4h")
    if (grepl("control group after 24 hours", tt)) return("control_24h")
  }
  if (block == "ES3") {
    if (grepl("20.*fa.*4 times", tt)) return("20uM_x4")
    if (grepl("50.*fa.*4 times", tt)) return("50uM_x4")
    if (grepl("control", tt)) return("control_x4")
  }
  NA_character_
}
md$condition <- mapply(classify_condition, md$title, md$block, USE.NAMES = FALSE)

if (nrow(md) != 39L || anyDuplicated(md$GSM) ||
    anyNA(md[, c("block", "source_ID", "bio_rep", "condition")])) {
  stop("Metadata audit failed: expected 39 uniquely mapped samples.")
}
if (!setequal(unique(md$block), c("ES1", "ES2", "ES3"))) {
  stop("Unexpected experiment blocks.")
}

# Confirm the GEO biological replicate label and source identifier are one-to-one
# within each block; every source/condition cell must contain exactly one array.
crosswalk <- unique(md[, c("block", "bio_rep", "source_ID")])
if (any(duplicated(crosswalk[, c("block", "bio_rep")])) ||
    any(duplicated(crosswalk[, c("block", "source_ID")]))) {
  stop("biological replicate and source ID do not map one-to-one within a block.")
}
coverage <- as.data.frame(table(md$block, md$source_ID, md$condition), stringsAsFactors = FALSE)
coverage <- coverage[coverage$Freq > 0, , drop = FALSE]
if (any(coverage$Freq != 1L)) stop("Duplicate or missing source-ID/condition observations.")

fwrite(md, file.path(OUT, "metadata", "GSE21477_metadata_sourceID_derived.csv"))
fwrite(crosswalk, file.path(OUT, "metadata", "GSE21477_sourceID_replicate_crosswalk.csv"))
fwrite(as.data.frame(table(md$block, md$condition)),
       file.path(OUT, "metadata", "GSE21477_block_condition_counts.csv"))

# ---------- 2. Download, unpack and map all 39 raw CEL files ----------
rawdir <- file.path(OUT, "rawCEL")
getGEOSuppFiles("GSE21477", makeDirectory = FALSE, baseDir = rawdir, fetch_files = TRUE)
raw_tar <- list.files(rawdir, pattern = "GSE21477_RAW\\.tar$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)
if (length(raw_tar) != 1L) {
  stop("Expected exactly one GSE21477_RAW.tar archive; inspect rawCEL directory.")
}
untar(raw_tar[1], exdir = rawdir)
cel_gz <- list.files(rawdir, pattern = "\\.CEL\\.gz$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
for (ff in cel_gz) {
  dest <- sub("\\.gz$", "", ff, ignore.case = TRUE)
  if (!file.exists(dest)) R.utils::gunzip(ff, destname = dest, remove = FALSE)
}
cel <- list.files(rawdir, pattern = "\\.CEL$", full.names = TRUE,
                  recursive = TRUE, ignore.case = TRUE)
gsm_from_file <- sub("^(GSM[0-9]+).*", "\\1", toupper(basename(cel)))
names(cel) <- gsm_from_file
if (length(cel) != 39L || anyDuplicated(gsm_from_file) || !all(md$GSM %in% names(cel))) {
  stop("CEL-to-GSM mapping failed: expected 39 unique CEL files matching all GSM accessions.")
}
md5 <- data.frame(GSM = names(cel), file = basename(cel), md5 = unname(tools::md5sum(cel)))
fwrite(md5, file.path(OUT, "metadata", "GSE21477_CEL_MD5.csv"))

# Keep array order identical to GEO metadata.
ab <- ReadAffy(filenames = unname(cel[md$GSM]))
sampleNames(ab) <- md$GSM
md <- md[match(sampleNames(ab), md$GSM), , drop = FALSE]
stopifnot(identical(sampleNames(ab), md$GSM))

# ---------- 3. Global RMA (common expression input for primary block models) ----------
es_global <- rma(ab)
expr_global <- exprs(es_global)
write.csv(data.frame(probe_id = rownames(expr_global), expr_global, check.names = FALSE),
          file.path(OUT, "global_RMA", "GSE21477_global_RMA_matrix.csv"), row.names = FALSE)

# ---------- 4. Contrasts ----------
contrast_pairs <- list(
  ES1_50uM_2h_vs_control_shared = c("50uM_2h", "control_shared"),
  ES1_50uM_4h_vs_control_shared = c("50uM_4h", "control_shared"),
  ES1_100uM_2h_vs_control_shared = c("100uM_2h", "control_shared"),
  ES2_100uM_4h_vs_control_4h = c("100uM_4h", "control_4h"),
  ES2_100uM_24h_vs_control_24h = c("100uM_24h", "control_24h"),
  ES2_200uM_4h_vs_control_4h = c("200uM_4h", "control_4h"),
  ES2_200uM_24h_vs_control_24h = c("200uM_24h", "control_24h"),
  ES3_20uM_x4_vs_control_x4 = c("20uM_x4", "control_x4"),
  ES3_50uM_x4_vs_control_x4 = c("50uM_x4", "control_x4")
)

# ---------- 5. Design-matrix helpers ----------
make_condition_matrix <- function(m) {
  conds <- sort(unique(as.character(m$condition)))
  z <- sapply(conds, function(cc) as.numeric(m$condition == cc))
  if (is.null(dim(z))) z <- matrix(z, ncol = 1L)
  colnames(z) <- paste0("cond_", make.names(conds))
  z
}

make_source_dummies_block <- function(m) {
  ids <- sort(unique(as.character(m$source_ID)))
  if (length(ids) < 2L) stop("At least two source IDs are required for blocking.")
  z <- sapply(ids[-1], function(id) as.numeric(m$source_ID == id))
  if (is.null(dim(z))) z <- matrix(z, ncol = 1L)
  colnames(z) <- paste0("src_", ids[-1])
  z
}

make_source_dummies_nested <- function(m) {
  pieces <- list()
  for (bb in sort(unique(m$block))) {
    ids <- sort(unique(as.character(m$source_ID[m$block == bb])))
    if (length(ids) < 2L) stop("At least two source IDs required in ", bb)
    for (id in ids[-1]) {
      nm <- paste0("src_", bb, "_", id)
      pieces[[nm]] <- as.numeric(m$block == bb & m$source_ID == id)
    }
  }
  do.call(cbind, pieces)
}

make_contrast_matrix <- function(design, pairs) {
  C <- matrix(0, nrow = ncol(design), ncol = length(pairs),
              dimnames = list(colnames(design), names(pairs)))
  for (j in seq_along(pairs)) {
    trt <- pairs[[j]][1]
    ctl <- pairs[[j]][2]
    trt_col <- paste0("cond_", make.names(trt))
    ctl_col <- paste0("cond_", make.names(ctl))
    if (!all(c(trt_col, ctl_col) %in% rownames(C))) {
      stop("Contrast condition missing from design: ", names(pairs)[j])
    }
    C[trt_col, j] <- 1
    C[ctl_col, j] <- -1
  }
  C
}

run_limma <- function(expr, m, pairs, prefix, model_label, design) {
  m <- m[match(colnames(expr), m$GSM), , drop = FALSE]
  if (!identical(colnames(expr), m$GSM)) stop("Expression and metadata order mismatch.")
  rank <- qr(design)$rank
  if (rank != ncol(design)) stop("Rank-deficient design for ", model_label)
  residual_df <- nrow(design) - rank
  if (residual_df < 2L) stop("Residual df <2 for ", model_label)
  C <- make_contrast_matrix(design, pairs)
  fit <- lmFit(expr, design)
  fit2 <- eBayes(contrasts.fit(fit, C), trend = TRUE, robust = TRUE)
  summary_list <- vector("list", ncol(C))
  for (j in seq_len(ncol(C))) {
    tt <- topTable(fit2, coef = j, number = Inf, adjust.method = "BH", sort.by = "none")
    tt$probe_id <- rownames(tt)
    tt$contrast <- colnames(C)[j]
    fwrite(tt, file.path(OUT, prefix, paste0(model_label, "_", make.names(colnames(C)[j]), ".csv")))
    summary_list[[j]] <- data.frame(
      model = model_label,
      contrast = colnames(C)[j],
      n_arrays = ncol(expr),
      design_rank = rank,
      residual_df = residual_df,
      FDR_sig_probe_sets = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
      FDR_sig_and_abs_logFC_ge_log2_1p5 = sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= log2(1.5), na.rm = TRUE),
      min_FDR = min(tt$adj.P.Val, na.rm = TRUE),
      median_abs_logFC = median(abs(tt$logFC), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  list(fit = fit2, summary = do.call(rbind, summary_list), design = design)
}

# ---------- 6. Primary analysis: separate block-specific models on common global RMA ----------
primary_summaries <- list()
primary_fits <- list()
for (bb in c("ES1", "ES2", "ES3")) {
  ii <- md$block == bb
  m_b <- md[ii, , drop = FALSE]
  expr_b <- expr_global[, m_b$GSM, drop = FALSE]
  pairs_b <- contrast_pairs[startsWith(names(contrast_pairs), paste0(bb, "_"))]
  design_b <- cbind(make_condition_matrix(m_b), make_source_dummies_block(m_b))
  out_b <- run_limma(expr_b, m_b, pairs_b, "global_RMA", paste0("PRIMARY_", bb), design_b)
  primary_summaries[[bb]] <- out_b$summary
  primary_fits[[bb]] <- out_b$fit
}
primary_summary <- do.call(rbind, primary_summaries)
fwrite(primary_summary, file.path(OUT, "global_RMA", "SUMMARY_PRIMARY_block_specific_sourceID_limma.csv"))

# ---------- 7. Sensitivity A: all-block source-ID nested fixed-effect model ----------
design_nested <- cbind(make_condition_matrix(md), make_source_dummies_nested(md))
all_pairs <- contrast_pairs
nested_out <- run_limma(expr_global, md, all_pairs, "global_RMA", "SENSITIVITY_global_sourceID_nested", design_nested)
fwrite(nested_out$summary,
       file.path(OUT, "global_RMA", "SUMMARY_SENSITIVITY_global_sourceID_nested.csv"))

# ---------- 8. Sensitivity B: historical all-block unblocked condition model ----------
design_unblocked <- make_condition_matrix(md)
unblocked_out <- run_limma(expr_global, md, all_pairs, "global_RMA", "SENSITIVITY_global_unblocked", design_unblocked)
fwrite(unblocked_out$summary,
       file.path(OUT, "global_RMA", "SUMMARY_SENSITIVITY_global_unblocked.csv"))

# ---------- 9. Sensitivity C: RMA separately within each experimental block ----------
within_summaries <- list()
for (bb in c("ES1", "ES2", "ES3")) {
  ii <- md$block == bb
  m_b <- md[ii, , drop = FALSE]
  ab_b <- ab[, m_b$GSM]
  es_b <- rma(ab_b)
  expr_b <- exprs(es_b)
  pairs_b <- contrast_pairs[startsWith(names(contrast_pairs), paste0(bb, "_"))]
  design_b <- cbind(make_condition_matrix(m_b), make_source_dummies_block(m_b))
  out_b <- run_limma(expr_b, m_b, pairs_b, "within_block_RMA", paste0("SENSITIVITY_within_RMA_", bb), design_b)
  within_summaries[[bb]] <- out_b$summary
  write.csv(data.frame(probe_id = rownames(expr_b), expr_b, check.names = FALSE),
            file.path(OUT, "within_block_RMA", paste0("GSE21477_", bb, "_RMA_matrix.csv")),
            row.names = FALSE)
}
within_summary <- do.call(rbind, within_summaries)
fwrite(within_summary,
       file.path(OUT, "within_block_RMA", "SUMMARY_SENSITIVITY_within_block_RMA.csv"))

# ---------- 10. Chip-level QC ----------
qc_written <- FALSE
if (requireNamespace("affyPLM", quietly = TRUE)) {
  suppressPackageStartupMessages(library(affyPLM))
  plm <- try(fitPLM(ab, background = TRUE, normalize = TRUE), silent = TRUE)
  if (!inherits(plm, "try-error")) {
    saveRDS(plm, file.path(OUT, "qc", "GSE21477_affyPLM.rds"))
    rle_mat <- RLE(plm)
    nuse_mat <- NUSE(plm)
    qc <- data.frame(
      GSM = sampleNames(ab),
      RLE_median = apply(rle_mat, 2, median, na.rm = TRUE),
      RLE_IQR = apply(rle_mat, 2, IQR, na.rm = TRUE),
      NUSE_median = apply(nuse_mat, 2, median, na.rm = TRUE),
      NUSE_IQR = apply(nuse_mat, 2, IQR, na.rm = TRUE)
    )
    qc <- merge(md, qc, by = "GSM", sort = FALSE)
    fwrite(qc, file.path(OUT, "qc", "GSE21477_RLE_NUSE_per_array.csv"))
    pdf(file.path(OUT, "qc", "GSE21477_RLE_NUSE_boxplots.pdf"), width = 11, height = 6)
    par(mfrow = c(1, 2))
    boxplot(rle_mat, las = 2, cex.axis = 0.55, main = "RLE by array", ylab = "RLE")
    abline(h = 0, lty = 2)
    boxplot(nuse_mat, las = 2, cex.axis = 0.55, main = "NUSE by array", ylab = "NUSE")
    abline(h = 1, lty = 2)
    dev.off()
    qc_written <- TRUE
  } else {
    writeLines(as.character(plm), file.path(OUT, "qc", "affyPLM_error.txt"))
  }
} else {
  writeLines("affyPLM unavailable; RLE/NUSE not generated.",
             file.path(OUT, "qc", "affyPLM_not_installed.txt"))
}

# PCA and expression distribution diagnostics from common global RMA.
vars <- matrixStats::rowVars(expr_global)
top <- order(vars, decreasing = TRUE)[seq_len(min(1000L, length(vars)))]
pc <- prcomp(t(expr_global[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca <- data.frame(GSM = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2])
pca <- merge(pca, md, by = "GSM", sort = FALSE)
fwrite(pca, file.path(OUT, "qc", "GSE21477_PCA_scores.csv"))
pdf(file.path(OUT, "qc", "GSE21477_PCA_by_block_and_condition.pdf"), width = 9, height = 7)
cols <- as.integer(factor(pca$block))
pch <- as.integer(factor(pca$condition)) + 14
plot(pca$PC1, pca$PC2, col = cols, pch = pch, xlab = "PC1", ylab = "PC2",
     main = "GSE21477 PCA: block (colour) and condition (symbol)")
text(pca$PC1, pca$PC2, labels = pca$GSM, pos = 3, cex = 0.5)
legend("topright", legend = levels(factor(pca$block)), col = seq_along(unique(pca$block)), pch = 19,
       title = "Block", cex = 0.75)
dev.off()

# Per-array expression distribution summaries.
expr_qc <- data.frame(
  GSM = colnames(expr_global),
  median_expression = apply(expr_global, 2, median, na.rm = TRUE),
  IQR_expression = apply(expr_global, 2, IQR, na.rm = TRUE),
  mean_expression = colMeans(expr_global, na.rm = TRUE)
)
expr_qc <- merge(md, expr_qc, by = "GSM", sort = FALSE)
fwrite(expr_qc, file.path(OUT, "qc", "GSE21477_expression_distribution_per_array.csv"))

# ---------- 11. Reproducibility record ----------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
writeLines(c(
  paste0("RLE_NUSE_generated=", qc_written),
  paste0("n_samples=", nrow(md)),
  paste0("n_probes=", nrow(expr_global)),
  paste0("primary_model=separate source-ID fixed-effect limma models within ES1/ES2/ES3"),
  paste0("sensitivity_models=global source-ID nested; global unblocked; within-block RMA"),
  "No manuscript result should be updated until all outputs and model assumptions are reviewed."
), file.path(OUT, "logs", "RUN_NOTES.txt"))

message("GSE21477 v2 reanalysis completed. Review all QC, design and summary outputs before manuscript revision.")
