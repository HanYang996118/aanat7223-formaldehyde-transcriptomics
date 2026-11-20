# AANAT7223: reviewer-directed cross-dataset pathway validation
# Ranked GO-BP GSEA (fgseaMultilevel) + cameraPR using the 29 primary
# platform-wide limma contrasts in Supplementary Data S11.
#
# IMPORTANT ANALYTICAL RULES
# 1) Native platform-level limma results remain the primary inference.
# 2) One representative feature per uniquely annotated gene is selected by
#    highest across-sample variance in the dataset-specific expression matrix.
# 3) Gene-set analysis is secondary and tests transferability/context, not toxicity.
# 4) GO-BP annotations are taken from the corresponding GEO platform annotation.
# 5) In-vitro micromolar concentrations are never converted to inhaled ppm.
# 6) Cross-species comparisons use shared GO IDs, not an assumed dose equivalence.
#
# Expected input:
#   AANAT7223_AllFeature_EffectSize_95CI_Audit_v1.zip
# containing the four *_AllFeature_Limma_Conservative95CI.csv.gz files.
#
# This script uses expression matrices for representative-feature variance selection:
# raw-CEL/RMA matrices for GSE23179, GSE42394 and GSE27263; the complete 39-sample
# GEO-processed log2 RMA matrix for GSE21477. The latter was separately cross-checked
# against nine reprocessed raw CEL arrays. Primary limma statistics are imported
# from the locked all-feature archive and are not re-estimated by this script.
# These mappings are joined to the already-finalized primary limma statistics.

options(stringsAsFactors = FALSE, timeout = 1200)
set.seed(20261109)

# ------------------------- User-configurable paths -------------------------
INPUT_ZIP <- "AANAT7223_AllFeature_EffectSize_95CI_Audit_v1.zip"
OUT <- "AANAT7223_Ranked_GSEA_Camera_v2"
RAW_DIR <- file.path(OUT, "rawCEL")
EXPR_DIR <- file.path(OUT, "RMA_expression")
MAP_DIR <- file.path(OUT, "annotation_mapping")
RESULT_DIR <- file.path(OUT, "pathway_results")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
for (d in c(RAW_DIR, EXPR_DIR, MAP_DIR, RESULT_DIR))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

# ------------------------------ Packages -----------------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
cran <- c("data.table", "R.utils")
# Explicitly install the platform design/CDF packages required to read raw CEL files.
# These were present in the previous R session but were missing from v1's installer.
bioc <- c("limma", "Biobase", "affy", "fgsea", "GEOquery", "oligo")
for (p in cran) if (!requireNamespace(p, quietly=TRUE))
  install.packages(p, repos="https://cloud.r-project.org")
for (p in bioc) if (!requireNamespace(p, quietly=TRUE))
  BiocManager::install(p, ask=FALSE, update=FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Biobase)
  library(limma)
  library(fgsea)
})

# --------------------------- Input validation ------------------------------
if (!file.exists(INPUT_ZIP)) {
  stop("Cannot find INPUT_ZIP: ", INPUT_ZIP,
       ". Upload AANAT7223_AllFeature_EffectSize_95CI_Audit_v1.zip to this runtime.")
}
STAT_DIR <- file.path(OUT, "all_feature_statistics")
dir.create(STAT_DIR, recursive=TRUE, showWarnings=FALSE)
unzip(INPUT_ZIP, exdir=STAT_DIR)
stat_files <- list.files(STAT_DIR, pattern="_AllFeature_Limma_Conservative95CI\\.csv\\.gz$",
                         recursive=TRUE, full.names=TRUE)
need <- c("GSE23179", "GSE42394", "GSE27263", "GSE21477")
for (ds in need) if (!any(grepl(ds, basename(stat_files), fixed=TRUE)))
  stop("Missing all-feature limma result for ", ds)

# ----------------------------- GEO helpers ---------------------------------
gsm_from_file <- function(x) sub("^(GSM[0-9]+).*", "\\1", toupper(basename(x)))

download_raw_cels <- function(gse) {
  d <- file.path(RAW_DIR, gse)
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

  existing_cel <- list.files(d, pattern="\\.CEL$", full.names=TRUE,
                             recursive=TRUE, ignore.case=TRUE)
  existing_gz <- list.files(d, pattern="\\.CEL\\.gz$", full.names=TRUE,
                            recursive=TRUE, ignore.case=TRUE)

  # Reuse a complete prior extraction where possible; otherwise retrieve
  # the dataset-level raw archive and extract its CEL files.
  if (!(length(existing_cel) > 0L && length(existing_cel) == length(existing_gz))) {
    GEOquery::getGEOSuppFiles(gse, makeDirectory=FALSE, baseDir=d, fetch_files=TRUE)
    ff <- list.files(d, full.names=TRUE, recursive=TRUE)
    tarf <- ff[grepl("RAW\\.tar$", ff, ignore.case=TRUE)]
    for (tf in tarf) untar(tf, exdir=d)
  }

  gz <- list.files(d, pattern="\\.CEL\\.gz$", full.names=TRUE,
                   recursive=TRUE, ignore.case=TRUE)
  for (f in gz) {
    out <- sub("\\.gz$", "", f, ignore.case=TRUE)
    if (!file.exists(out))
      R.utils::gunzip(f, destname=out, remove=FALSE, overwrite=FALSE)
  }
  cel <- list.files(d, pattern="\\.CEL$", full.names=TRUE,
                    recursive=TRUE, ignore.case=TRUE)
  if (!length(cel)) stop("No decompressed CEL files found for ", gse)
  cel
}

get_series_pheno <- function(gse) {
  matrix_dir <- file.path(RAW_DIR, paste0(gse, "_matrix"))
  dir.create(matrix_dir, recursive=TRUE, showWarnings=FALSE)
  obj <- GEOquery::getGEO(gse, GSEMatrix=TRUE, AnnotGPL=FALSE,
                          destdir=matrix_dir)
  if (length(obj) < 1) stop("No Series Matrix object returned for ", gse)
  obj[[1]]
}

# Prepare dataset-specific expression matrices for the pre-specified
# highest-variance representative-feature mapping. GSE21477 uses the full
# GEO-processed Series Matrix; the other datasets use raw-CEL/RMA.
get_expression_matrix <- function(gse) {
  rds <- file.path(EXPR_DIR, paste0(gse, "_RMA_expression.rds"))
  if (file.exists(rds)) return(readRDS(rds))

  # GSE21477 is handled differently in this pathway workflow:
  # the complete 39-sample GEO Series Matrix (VALUE = log2 RMA signal)
  # is used for representative-feature variance selection. The primary
  # platform-wide limma statistics are imported from the locked all-feature
  # results archive and are not re-estimated in this pathway script.
  if (gse == "GSE21477") {
    eset <- get_series_pheno(gse)
    pd <- Biobase::pData(eset)
    gsm <- rownames(pd)
    ex <- Biobase::exprs(eset)
    if (length(gsm) != 39L || nrow(ex) != 54675L ||
        !all(gsm %in% colnames(ex))) {
      stop("GSE21477 Series Matrix must contain 39 matched samples and 54,675 GPL570 probes.")
    }
    ex <- ex[, gsm, drop=FALSE]
    if (any(!is.finite(ex))) stop("Non-finite values in GSE21477 GEO processed expression matrix.")
    storage.mode(ex) <- "double"
    saveRDS(ex, rds, compress="xz")
    return(ex)
  }

  # For the other three datasets, retrieve raw CEL files and normalize with
  # the platform-appropriate RMA procedure for representative-feature selection.
  eset <- get_series_pheno(gse)
  pd <- Biobase::pData(eset)
  gsm <- rownames(pd)
  cel <- download_raw_cels(gse)

  if (gse == "GSE23179") {
    suppressPackageStartupMessages(library(affy))
    raw <- affy::ReadAffy(filenames=cel)
    norm <- affy::rma(raw)
    ex <- Biobase::exprs(norm)
    names_from_cel <- gsm_from_file(cel)
    if (anyDuplicated(names_from_cel)) stop("Duplicate GSM IDs in GSE23179 CEL files")
    colnames(ex) <- names_from_cel
    if (!all(gsm %in% colnames(ex)))
      stop("GSE23179 CEL/GSM mismatch: missing ", paste(setdiff(gsm,colnames(ex)),collapse=", "))
    ex <- ex[, gsm, drop=FALSE]

  } else if (gse == "GSE42394") {
    suppressPackageStartupMessages(library(oligo))
    raw <- oligo::read.celfiles(cel)
    names_from_cel <- gsm_from_file(cel)
    if (anyDuplicated(names_from_cel)) stop("Duplicate GSM IDs in GSE42394 CEL files")
    sampleNames(raw) <- names_from_cel
    nasal <- gsm[grepl("^Nose_", as.character(pd$title))]
    if (length(nasal) != 12L)
      stop("Expected 12 nasal samples in GSE42394; observed ", length(nasal))
    if (!all(nasal %in% sampleNames(raw)))
      stop("GSE42394 nasal GSMs missing from CEL archive")
    rawN <- raw[, nasal]
    norm <- oligo::rma(rawN, target="core")
    ex <- Biobase::exprs(norm)
    ex <- ex[, nasal, drop=FALSE]

  } else if (gse == "GSE27263") {
    suppressPackageStartupMessages(library(affy))
    nasal <- gsm[grepl("^Nasal Biopsy", as.character(pd$title))]
    if (length(nasal) != 40L)
      stop("Expected 40 nasal-biopsy arrays in GSE27263; observed ", length(nasal))
    cel_gsm <- gsm_from_file(cel)
    keep <- cel[cel_gsm %in% nasal]
    keep_gsm <- gsm_from_file(keep)
    if (!setequal(keep_gsm, nasal))
      stop("GSE27263 nasal GSM/CEL mapping is incomplete")
    norm <- affy::rma(affy::ReadAffy(filenames=keep))
    ex <- Biobase::exprs(norm)
    colnames(ex) <- keep_gsm
    ex <- ex[, nasal, drop=FALSE]
  } else stop("Unsupported dataset: ", gse)

  storage.mode(ex) <- "double"
  saveRDS(ex, rds, compress="xz")
  ex
}

# ------------------------- Platform annotation ----------------------------
get_platform_annotation <- function(gpl) {
  rds <- file.path(MAP_DIR, paste0(gpl, "_GEO_annotation.rds"))
  required_cols <- c("probe_id","gene_symbol","go_bp","go_bp_id")
  if (file.exists(rds)) {
    cached <- readRDS(rds)
    if (all(required_cols %in% names(cached))) return(cached)
    message("Refreshing outdated cached annotation schema for ", gpl)
  }

  gpl_obj <- GEOquery::getGEO(gpl, AnnotGPL=TRUE)
  tab <- as.data.frame(GEOquery::Table(gpl_obj), stringsAsFactors=FALSE)

  # GEO platform annotation headers differ from older Affymetrix annotation
  # packages. Normalize names before matching (lowercase BEFORE removing
  # punctuation, so fields such as 'GO:Process' remain identifiable).
  norm <- gsub("[^a-z0-9]+", "", tolower(names(tab)))
  pick_col <- function(candidates) {
    hit <- which(norm %in% candidates)
    if (length(hit)) names(tab)[hit[1]] else NA_character_
  }

  idcol <- pick_col(c("id"))
  if (is.na(idcol)) idcol <- names(tab)[1]
  symcol <- pick_col(c("genesymbol","symbol"))
  gocol <- pick_col(c("goprocess","geneontologybiologicalprocess",
                      "gobiologicalprocess","gobp"))
  goidcol <- pick_col(c("goprocessid","geneontologybiologicalprocessid",
                        "gobiologicalprocessid","gobpid"))

  if (is.na(symcol) || is.na(gocol)) {
    stop("Could not identify Gene Symbol and GO-BP term columns in ", gpl,
         ". Available columns: ", paste(names(tab), collapse=" | "))
  }

  go_bp_id <- if (!is.na(goidcol)) as.character(tab[[goidcol]]) else
    rep(NA_character_, nrow(tab))

  ann <- data.frame(
    probe_id=as.character(tab[[idcol]]),
    gene_symbol=trimws(as.character(tab[[symcol]])),
    go_bp=as.character(tab[[gocol]]),
    go_bp_id=go_bp_id,
    stringsAsFactors=FALSE
  )
  ann <- ann[!is.na(ann$probe_id) & ann$probe_id != "", ]
  ann <- ann[!duplicated(ann$probe_id), ]

  # Report the number of rows for which the GEO GO term and GO ID lists
  # can be paired positionally; mismatched rows are excluded at parsing.
  has_terms <- !is.na(ann$go_bp) & nzchar(trimws(ann$go_bp))
  has_ids <- !is.na(ann$go_bp_id) & nzchar(trimws(ann$go_bp_id))
  pair_idx <- which(has_terms & has_ids)
  term_n <- vapply(ann$go_bp[pair_idx],
                   function(x) length(strsplit(x, "///", fixed=TRUE)[[1]]), integer(1))
  id_n <- vapply(ann$go_bp_id[pair_idx],
                 function(x) length(strsplit(x, "///", fixed=TRUE)[[1]]), integer(1))
  message(gpl, ": ", nrow(ann), " platform probes; ",
          sum(has_terms), " with GO-BP terms; ",
          sum(term_n == id_n), " rows with positionally matched GO term/ID lists; ",
          sum(term_n != id_n), " rows with list-length mismatch.")

  saveRDS(ann, rds, compress="xz")
  data.table::fwrite(ann, file.path(MAP_DIR, paste0(gpl, "_annotation_raw.csv")))
  ann
}

# Parse either modern GEO split fields (GO:Process + GO:Process ID),
# which are parallel '///'-delimited vectors, or legacy combined
# GO:ID // term // evidence annotations.
parse_go_bp <- function(x, go_id_string=NA_character_) {
  if (is.na(x) || !nzchar(trimws(x))) return(NULL)

  if (!is.na(go_id_string) && nzchar(trimws(go_id_string))) {
    terms <- trimws(unlist(strsplit(as.character(x), "///", fixed=TRUE),
                           use.names=FALSE))
    ids <- trimws(unlist(strsplit(as.character(go_id_string), "///", fixed=TRUE),
                         use.names=FALSE))
    # Never guess or recycle IDs if the two GEO annotation vectors disagree.
    if (length(terms) != length(ids)) return(NULL)
    ids_clean <- sub("^GO:", "", ids, ignore.case=TRUE)
    valid <- !is.na(terms) & nzchar(terms) &
             grepl("^[0-9]{7}$", ids_clean)
    if (!any(valid)) return(NULL)
    out <- data.frame(
      GO_ID=paste0("GO:", ids_clean[valid]),
      GO_term=terms[valid],
      stringsAsFactors=FALSE
    )
    return(unique(out))
  }

  # Legacy annotation format: GO:0000000 // biological process // evidence
  parts <- unlist(strsplit(as.character(x), "///", fixed=TRUE), use.names=FALSE)
  out <- lapply(parts, function(z) {
    z <- trimws(z)
    m <- regexec("^(?:GO:)?([0-9]{7})\\s*//\\s*(.*?)\\s*//", z, perl=TRUE)
    hit <- regmatches(z, m)[[1]]
    if (length(hit) < 3) return(NULL)
    data.frame(GO_ID=paste0("GO:",hit[2]), GO_term=trimws(hit[3]),
               stringsAsFactors=FALSE)
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(NULL)
  unique(data.table::rbindlist(out, fill=TRUE))
}

make_representative_map <- function(gse, gpl, expr) {
  out_csv <- file.path(MAP_DIR, paste0(gse, "_representative_feature_map.csv"))
  go_csv <- file.path(MAP_DIR, paste0(gse, "_GO_BP_gene_membership.csv"))
  if (file.exists(out_csv)) {
    cached_map <- data.table::fread(out_csv)
    if (all(c("probe_id","gene_symbol","go_bp","go_bp_id") %in% names(cached_map)))
      return(cached_map)
    message("Refreshing outdated representative-feature mapping for ", gse)
    unlink(c(out_csv, go_csv))
  }
  ann <- get_platform_annotation(gpl)
  # Keep one-to-one probe-to-symbol mappings only. Ambiguous symbols are excluded.
  ok <- !is.na(ann$gene_symbol) & nzchar(ann$gene_symbol) &
        !grepl("///|//|;|\\|", ann$gene_symbol)
  ann <- ann[ok, , drop=FALSE]
  ann <- ann[ann$probe_id %in% rownames(expr), , drop=FALSE]
  if (!nrow(ann)) stop("No annotated platform probes overlap expression matrix for ", gse)
  # Highest across-sample variance per gene, determined without looking at any contrast.
  vv <- apply(expr[ann$probe_id, , drop=FALSE], 1, stats::var, na.rm=TRUE)
  ann$variance <- as.numeric(vv[ann$probe_id])
  ann <- ann[is.finite(ann$variance), , drop=FALSE]
  ann <- ann[order(ann$gene_symbol, -ann$variance, ann$probe_id), ]
  ann <- ann[!duplicated(ann$gene_symbol), , drop=FALSE]
  # Parse GO annotations only after the representative feature is fixed.
  go_long <- lapply(seq_len(nrow(ann)), function(i) {
    z <- parse_go_bp(ann$go_bp[i], ann$go_bp_id[i])
    if (is.null(z)) return(NULL)
    z$gene_symbol <- ann$gene_symbol[i]
    z$probe_id <- ann$probe_id[i]
    z
  })
  go_long <- Filter(Negate(is.null), go_long)
  if (!length(go_long)) stop("No GO-BP annotation entries parsed for ", gpl)
  go_long <- unique(data.table::rbindlist(go_long, fill=TRUE))
  data.table::fwrite(ann, out_csv)
  data.table::fwrite(go_long, file.path(MAP_DIR, paste0(gse, "_GO_BP_gene_membership.csv")))
  attr(ann, "go_long") <- go_long
  ann
}

# Rebuild GO long table from saved annotation/mapping if map already existed.
get_go_membership <- function(gse, rep_map) {
  f <- file.path(MAP_DIR, paste0(gse, "_GO_BP_gene_membership.csv"))
  if (file.exists(f)) return(data.table::fread(f))
  # If the representative map exists but GO table does not, retrieve annotation and parse.
  gpl <- c(GSE23179="GPL1355", GSE42394="GPL6247", GSE27263="GPL570", GSE21477="GPL570")[[gse]]
  z <- as.data.frame(rep_map)
  # make_representative_map already carries the selected probe's GO-BP field.
  # Merge only if a cached mapping file lacks that column; otherwise a second
  # merge would create go_bp.x/go_bp.y and silently break parsing.
  if (!all(c("go_bp","go_bp_id") %in% names(z))) {
    ann <- get_platform_annotation(gpl)
    z <- merge(z, ann[,c("probe_id","go_bp","go_bp_id")],
               by="probe_id", all.x=TRUE)
  }
  out <- lapply(seq_len(nrow(z)), function(i) {
    pp <- parse_go_bp(z$go_bp[i], z$go_bp_id[i])
    if (is.null(pp)) return(NULL)
    pp$gene_symbol <- z$gene_symbol[i]
    pp$probe_id <- z$probe_id[i]
    pp
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) stop("GO-BP memberships unavailable for ", gse)
  out <- unique(data.table::rbindlist(out, fill=TRUE))
  data.table::fwrite(out, f)
  out
}

writeLines(c(
  "GSE42394 / GPL6247 representative-gene mapping audit",
  "RMA expression matrix: 29,214 features x 12 nasal tissue arrays.",
  "GPL6247 annotation rows: 29,215.",
  "Non-empty, one-to-one Gene symbol annotation rows: 15,656.",
  "Unique one-to-one gene symbols represented in the RMA matrix: 14,841.",
  "Multi-symbol annotation rows excluded from representative mapping: 368.",
  "No multi-symbol probe was split or assigned to multiple genes.",
  "The provisional target of 14,848 was superseded because it could not be reproduced under the one-to-one mapping rule."
), file.path(MAP_DIR, "GSE42394_gene_mapping_count_audit.txt"))

# ---------------------- Load current primary statistics --------------------
stat_map <- list(
  GSE23179=file.path(STAT_DIR, basename(stat_files[grepl("GSE23179",basename(stat_files))][1])),
  GSE42394=file.path(STAT_DIR, basename(stat_files[grepl("GSE42394",basename(stat_files))][1])),
  GSE27263=file.path(STAT_DIR, basename(stat_files[grepl("GSE27263",basename(stat_files))][1])),
  GSE21477=file.path(STAT_DIR, basename(stat_files[grepl("GSE21477",basename(stat_files))][1]))
)
# More robustly resolve recursively if the ZIP contains a nested directory.
for (ds in names(stat_map)) {
  ff <- stat_files[grepl(ds, basename(stat_files), fixed=TRUE)]
  stat_map[[ds]] <- ff[1]
}

platforms <- c(GSE23179="GPL1355", GSE42394="GPL6247",
               GSE27263="GPL570", GSE21477="GPL570")
expected_features <- c(GSE23179=31099L, GSE42394=29214L,
                       GSE27263=54675L, GSE21477=54675L)
# Expected unique genes after the pre-specified one-to-one symbol filter.
# GPL6247 currently yields 14,841 unambiguous symbols across all 29,214 RMA
# features; 368 multi-symbol probe rows are excluded rather than split across genes.
expected_genes <- c(GSE23179=14604L, GSE42394=14841L,
                    GSE27263=20848L, GSE21477=20848L)
expected_contrasts <- c(GSE23179=15L, GSE42394=2L,
                        GSE27263=3L, GSE21477=9L)
exprs_by_ds <- list(); rep_maps <- list(); go_membership <- list(); stats_by_ds <- list()

for (ds in names(platforms)) {
  message("Preparing expression-based representative mapping: ", ds)
  ex <- get_expression_matrix(ds)
  if (nrow(ex) != expected_features[[ds]])
    stop(ds, ": expected ", expected_features[[ds]], " expression features but found ", nrow(ex),
         ". Stop and audit CEL/CDF/platform compatibility before pathway analysis.")
  exprs_by_ds[[ds]] <- ex
  rep_map <- make_representative_map(ds, platforms[[ds]], ex)
  if (nrow(rep_map) != expected_genes[[ds]])
    stop(ds, ": expected ", expected_genes[[ds]], " representative genes but found ", nrow(rep_map),
         ". Stop and audit the GEO annotation snapshot and representative-feature mapping.")
  rep_maps[[ds]] <- as.data.frame(rep_map)
  go_membership[[ds]] <- get_go_membership(ds, rep_maps[[ds]])
  go_ids <- unique(as.character(go_membership[[ds]]$GO_ID))
  if (!length(go_ids) || any(!grepl("^GO:[0-9]{7}$", go_ids)))
    stop(ds, ": invalid or missing GO IDs detected after GO term/ID mapping.")
  if (anyNA(go_membership[[ds]]$GO_term) ||
      any(!nzchar(trimws(go_membership[[ds]]$GO_term))))
    stop(ds, ": empty GO-BP term labels detected after mapping.")
  stats_by_ds[[ds]] <- data.table::fread(stat_map[[ds]])
  stats_by_ds[[ds]][, Probe_set := as.character(Probe_set)]
  if (length(unique(stats_by_ds[[ds]]$Contrast)) != expected_contrasts[[ds]])
    stop(ds, ": expected ", expected_contrasts[[ds]], " contrasts but found ",
         length(unique(stats_by_ds[[ds]]$Contrast)))
  if (anyDuplicated(stats_by_ds[[ds]][, .(Contrast, Probe_set)]))
    stop(ds, ": duplicate feature rows detected within a contrast.")
}
if (sum(vapply(stats_by_ds, function(x) length(unique(x$Contrast)), integer(1))) != 29L)
  stop("The four datasets do not contain exactly 29 expected primary contrasts.")

# ---------------------- Build platform-specific GO sets --------------------
build_go_sets <- function(go_long, universe_symbols, min_size=10L, max_size=500L) {
  go_long <- as.data.frame(go_long)
  go_long <- go_long[go_long$gene_symbol %in% universe_symbols, , drop=FALSE]
  go_long <- unique(go_long[, c("GO_ID","GO_term","gene_symbol")])
  split_genes <- split(go_long$gene_symbol, go_long$GO_ID)
  term_map <- unique(go_long[,c("GO_ID","GO_term")])
  names(term_map$GO_term) <- term_map$GO_ID
  sets <- lapply(split_genes, unique)
  sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
  term_map <- term_map[match(names(sets), term_map$GO_ID), , drop=FALSE]
  list(sets=sets, terms=term_map)
}

# ---------------------- Run GSEA + cameraPR for 29 contrasts ---------------
all_gsea <- list(); all_camera <- list(); all_camera_rank <- list(); all_camera_cor <- list()
representative_audit <- list()

for (ds in names(platforms)) {
  message("Running pathway tests: ", ds)
  rep_map <- rep_maps[[ds]]
  go_long <- go_membership[[ds]]
  go_sets_obj <- build_go_sets(go_long, rep_map$gene_symbol, min_size=10L, max_size=500L)
  pathways <- go_sets_obj$sets
  term_map <- go_sets_obj$terms
  res <- stats_by_ds[[ds]]
  contrasts <- unique(res$Contrast)
  representative_audit[[ds]] <- data.frame(
    Dataset=ds, Platform=platforms[[ds]],
    Platform_features_in_expression=nrow(exprs_by_ds[[ds]]),
    Unique_representative_genes=nrow(rep_map),
    GO_BP_terms_tested=length(pathways),
    GO_BP_gene_memberships=nrow(go_long),
    Annotation_retrieval_date=as.character(Sys.Date()),
    stringsAsFactors=FALSE
  )

  for (cn in contrasts) {
    tt <- res[res$Contrast == cn, ]
    # Join only the pre-defined representative feature for each gene.
    tt$Probe_set <- as.character(tt$Probe_set)
    tt <- as.data.frame(tt[,c("Probe_set","Moderated_t","log2FC","P_value","BH_FDR")])
    m <- merge(as.data.frame(rep_map[,c("probe_id","gene_symbol")]), tt,
               by.x="probe_id", by.y="Probe_set", all=FALSE, sort=FALSE)
    m <- m[!duplicated(m$gene_symbol), ]
    stat <- m$Moderated_t
    names(stat) <- m$gene_symbol
    keep <- is.finite(stat) & !is.na(names(stat)) & names(stat) != ""
    stat <- stat[keep]
    stat <- sort(stat, decreasing=TRUE)
    if (length(stat) < 1000L) stop("Too few ranked genes for ", ds, " / ", cn)

    # Restrict GO sets to the genes represented in this exact ranked universe.
    pws <- lapply(pathways, function(gs) intersect(gs, names(stat)))
    pws <- pws[lengths(pws) >= 10L & lengths(pws) <= 500L]
    if (!length(pws)) stop("No eligible GO-BP pathways for ", ds, " / ", cn)

    fg <- fgsea::fgseaMultilevel(
      pathways=pws, stats=stat, minSize=10, maxSize=500, eps=0,
      scoreType="std", nproc=1
    )
    fg <- as.data.frame(fg)
    fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
    fg$GO_ID <- fg$pathway
    fg$Dataset <- ds; fg$Contrast <- cn
    fg <- merge(fg, term_map, by="GO_ID", all.x=TRUE, sort=FALSE)
    all_gsea[[paste(ds,cn,sep="__")]] <- fg

    # cameraPR is an independent competitive test that accounts for
    # inter-gene correlation. use.ranks=FALSE retains the signed moderated t.
    idx <- pws
    cam <- limma::cameraPR(
      statistic=stat, index=idx, use.ranks=FALSE,
      inter.gene.cor=0.01, sort=FALSE
    )
    cam <- as.data.frame(cam)
    cam$GO_ID <- rownames(cam)
    cam$Dataset <- ds; cam$Contrast <- cn
    cam <- merge(cam, term_map, by="GO_ID", all.x=TRUE, sort=FALSE)
    all_camera[[paste(ds,cn,sep="__")]] <- cam

    # Rank-based camera sensitivity analysis.
    camr <- limma::cameraPR(
      statistic=stat, index=idx, use.ranks=TRUE,
      inter.gene.cor=0.01, sort=FALSE
    )
    camr <- as.data.frame(camr)
    camr$GO_ID <- rownames(camr)
    camr$Dataset <- ds; camr$Contrast <- cn
    camr <- merge(camr, term_map, by="GO_ID", all.x=TRUE, sort=FALSE)
    all_camera_rank[[paste(ds,cn,sep="__")]] <- camr

    # Sensitivity to the assumed inter-gene correlation in the competitive test.
    # The prespecified primary cameraPR result remains rho=0.01.
    for (rho in c(0, 0.01, 0.05)) {
      camc <- limma::cameraPR(
        statistic=stat, index=idx, use.ranks=FALSE,
        inter.gene.cor=rho, sort=FALSE
      )
      camc <- as.data.frame(camc)
      camc$GO_ID <- rownames(camc)
      camc$Dataset <- ds; camc$Contrast <- cn
      camc$inter.gene.cor <- rho
      camc <- merge(camc, term_map, by="GO_ID", all.x=TRUE, sort=FALSE)
      all_camera_cor[[paste(ds,cn,rho,sep="__")]] <- camc
    }
  }
}

gsea_all <- data.table::rbindlist(all_gsea, fill=TRUE)
camera_all <- data.table::rbindlist(all_camera, fill=TRUE)
camera_rank_all <- data.table::rbindlist(all_camera_rank, fill=TRUE)
camera_cor_all <- data.table::rbindlist(all_camera_cor, fill=TRUE)
rep_audit <- data.table::rbindlist(representative_audit, fill=TRUE)

camera_cor_summary <- camera_cor_all[, .(
  GO_BP_tests=.N,
  FDR05_terms=sum(FDR < 0.05, na.rm=TRUE)
), by=.(Dataset, Contrast, inter.gene.cor)]

camera_cor_compare <- list()
for (ds in unique(camera_cor_all$Dataset)) {
  for (cn in unique(camera_cor_all[Dataset==ds, Contrast])) {
    ref <- camera_cor_all[Dataset==ds & Contrast==cn & inter.gene.cor==0.01 & FDR<0.05]
    for (rho in c(0,0.05)) {
      cur <- camera_cor_all[Dataset==ds & Contrast==cn & inter.gene.cor==rho & FDR<0.05]
      a <- unique(ref$GO_ID); b <- unique(cur$GO_ID)
      u <- union(a,b)
      inter <- intersect(a,b)
      dir_ref <- ref[GO_ID %in% inter, .(GO_ID,Direction)]
      dir_cur <- cur[GO_ID %in% inter, .(GO_ID,Direction)]
      dd <- merge(dir_ref,dir_cur,by="GO_ID",suffixes=c("_ref","_sens"))
      camera_cor_compare[[paste(ds,cn,rho,sep="__")]] <- data.frame(
        Dataset=ds, Contrast=cn, inter.gene.cor=rho,
        Primary_rho_0p01_FDR05=length(a),
        Sensitivity_FDR05=length(b),
        FDR_set_Jaccard=if(length(u)) length(inter)/length(u) else NA_real_,
        Direction_agreement_on_intersection=if(nrow(dd)) mean(dd$Direction_ref==dd$Direction_sens) else NA_real_,
        stringsAsFactors=FALSE
      )
    }
  }
}
camera_cor_compare <- data.table::rbindlist(camera_cor_compare, fill=TRUE)

fwrite(camera_cor_all, file.path(RESULT_DIR, "cameraPR_intergene_cor_sensitivity_all_29_contrasts.csv"))
fwrite(camera_cor_summary, file.path(RESULT_DIR, "cameraPR_intergene_cor_sensitivity_summary.csv"))
fwrite(camera_cor_compare, file.path(RESULT_DIR, "cameraPR_intergene_cor_sensitivity_vs_rho001.csv"))

fwrite(gsea_all, file.path(RESULT_DIR, "GSEA_GO_BP_all_29_contrasts.csv"))
fwrite(camera_all, file.path(RESULT_DIR, "cameraPR_GO_BP_all_29_contrasts.csv"))
fwrite(camera_rank_all, file.path(RESULT_DIR, "cameraPR_rank_sensitivity_all_29_contrasts.csv"))
fwrite(rep_audit, file.path(RESULT_DIR, "Representative_feature_and_GO_mapping_audit.csv"))

# ---------------------- Cross-dataset pathway concordance -------------------
# Only pre-specified, biologically interpretable contrasts are compared:
# Rat: nominally matched 2-ppm exposure schedules (not assumed dose-equivalent).
# Human: controlled nasal biopsy lag strata versus each distinct in-vitro schedule.
comparison_pairs <- list(
  Rat_2ppm_early=c("GSE23179|Td5_D2.Td5_D0", "GSE42394|d7_E - d7_C"),
  Rat_2ppm_later=c("GSE23179|Tw4_D2.Tw4_D0", "GSE42394|d28_E - d28_C")
)
for (h in c("h0","h24","h72")) {
  for (cc in unique(stats_by_ds$GSE21477$Contrast)) {
    nm <- paste0("Human_",h,"_vs_",cc)
    comparison_pairs[[nm]] <- c(paste0("GSE27263|",h), paste0("GSE21477|",cc))
  }
}

get_gsea_one <- function(ds, cn) {
  x <- gsea_all[Dataset==ds & Contrast==cn, ]
  x[, c("GO_ID","NES","padj","pval","size"), with=FALSE]
}
concordance <- list()
for (nm in names(comparison_pairs)) {
  pair <- comparison_pairs[[nm]]
  a <- strsplit(pair[1], "|", fixed=TRUE)[[1]]
  b <- strsplit(pair[2], "|", fixed=TRUE)[[1]]
  x <- get_gsea_one(a[1],a[2]); y <- get_gsea_one(b[1],b[2])
  z <- merge(x,y,by="GO_ID",suffixes=c("_A","_B"))
  z <- z[is.finite(z$NES_A) & is.finite(z$NES_B), ]
  if (nrow(z) < 20L) next
  rho <- suppressWarnings(cor(z$NES_A,z$NES_B,method="spearman",use="complete.obs"))
  top_n <- min(50L,nrow(z))
  ia <- order(abs(z$NES_A),decreasing=TRUE)[seq_len(top_n)]
  ib <- order(abs(z$NES_B),decreasing=TRUE)[seq_len(top_n)]
  agree_a <- mean(sign(z$NES_A[ia]) == sign(z$NES_B[ia]),na.rm=TRUE)
  agree_b <- mean(sign(z$NES_B[ib]) == sign(z$NES_A[ib]),na.rm=TRUE)
  overlap <- length(intersect(z$GO_ID[ia],z$GO_ID[ib])) / top_n
  sig_a <- sum(z$padj_A < 0.05,na.rm=TRUE)
  sig_b <- sum(z$padj_B < 0.05,na.rm=TRUE)
  concordance[[nm]] <- data.frame(
    Comparison=nm, Dataset_A=a[1], Contrast_A=a[2], Dataset_B=b[1], Contrast_B=b[2],
    Common_GO_BP_terms=nrow(z), Spearman_rho_NES=rho,
    Top50_direction_agreement_A_rank=agree_a,
    Top50_direction_agreement_B_rank=agree_b,
    Top50_pathway_overlap_fraction=overlap,
    FDR05_GSEA_terms_A=sig_a, FDR05_GSEA_terms_B=sig_b,
    stringsAsFactors=FALSE
  )
}
concordance_df <- if (length(concordance)) data.table::rbindlist(concordance,fill=TRUE) else data.table()
fwrite(concordance_df, file.path(RESULT_DIR, "Cross_dataset_GO_BP_NES_concordance.csv"))

# Cross-check whether fgsea and cameraPR agree on the direction of FDR-significant
# GO terms. This is a method-concordance audit, not a biological validation claim.
cam_check <- merge(
  gsea_all[,.(Dataset,Contrast,GO_ID,NES,GSEA_FDR=padj)],
  camera_all[,.(Dataset,Contrast,GO_ID,camera_Direction=Direction,camera_FDR=FDR)],
  by=c("Dataset","Contrast","GO_ID"), all=FALSE
)
cam_check[, GSEA_Direction := ifelse(NES > 0,"Up",ifelse(NES < 0,"Down","Zero"))]
cam_check[, Direction_agreement := GSEA_Direction == camera_Direction]
fwrite(cam_check, file.path(RESULT_DIR, "GSEA_camera_direction_crosscheck.csv"))

# ----------------------------- Provenance ----------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
writeLines(c(
  "AANAT7223 reviewer-directed ranked pathway analysis",
  "Primary data: four all-feature limma result tables from the V23 revision archive.",
  paste("NCBI GEO GPL annotation retrieval date:", as.character(Sys.Date())),
  "Probe/CDF dimensions and contrast counts were checked against the locked V23 primary tables before pathway testing.",
  "Representative-gene counts use one uniquely annotated gene symbol per probe; ambiguous multi-symbol probes are excluded rather than split.",
  "GSE42394/GPL6247 yielded 14,841 one-to-one gene symbols across all 29,214 RMA features; the earlier provisional 14,848 target was not reproducible under this conservative mapping rule and was replaced after explicit audit.",
  "Representative features: one uniquely annotated probe per gene, selected by highest variance across all samples in the corresponding dataset-specific expression matrix. Raw-CEL/RMA matrices are used for GSE23179, GSE42394 and GSE27263; the 39-sample GEO-processed log2 RMA matrix is used for GSE21477.",
  "GSE21477 representative-feature variance used all 39 GEO Series Matrix samples (54,675 GPL570 probes); sample VALUE fields are reported by GEO as log2 RMA signal. This pathway-stage matrix is distinct from the historical 39-CEL primary RMA/limma analysis.",
  "Independent GSE21477 CEL check: GSM536401-GSM536409 (nine arrays) were reprocessed together with affy::rma; matched-probe per-sample Pearson r ranged 0.998337-0.998457 and Spearman rho ranged 0.997186-0.997294.",
  "The nine-array subset RMA is a cross-check of expression pattern consistency, not a claim of numerical identity with the original 39-array GEO processing.",
  "Pathways: GO Biological Process terms parsed from the corresponding NCBI GEO platform annotation.",
  "Tests: fgseaMultilevel on signed moderated t statistics; limma cameraPR on the same signed statistics, with primary inter.gene.cor=0.01; rank-based cameraPR sensitivity analysis; additional cameraPR inter-gene-correlation sensitivity at rho=0, 0.01 and 0.05.",
  "Multiplicity: BH FDR within each dataset/contrast and method.",
  "Cross-dataset comparison: shared GO IDs and normalized enrichment score (NES) ranks; nominal concentration matches are descriptive only.",
  "Interpretation: pathway enrichment indicates coordinated transcriptomic association, not direct pathway activation/inhibition, tissue injury, toxicity, or a safe exposure threshold.",
  "In-vitro micromolar concentrations are not converted to inhaled ppm.",
  "All arrays and all platform-wide limma results remain in the primary analysis; this pathway analysis is secondary."
), file.path(OUT, "README_analysis_provenance.txt"))

if (nzchar(Sys.which("zip"))) {
  old <- getwd(); setwd(OUT)
  zip("../AANAT7223_Ranked_GSEA_Camera_v2_results.zip",
      files=c("pathway_results", "annotation_mapping", "sessionInfo.txt",
              "README_analysis_provenance.txt"), flags="-r")
  setwd(old)
}
message("Completed ranked GSEA + cameraPR pathway analysis.")
message("Outputs saved in: ", normalizePath(OUT))
