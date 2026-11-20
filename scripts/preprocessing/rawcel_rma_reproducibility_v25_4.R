#!/usr/bin/env Rscript
# AANAT7223 — raw-CEL/RMA reproducibility for the three GEO datasets used
# for representative-feature variance selection.
# This script does not refit the primary all-feature limma models.
# It uses the locally archived GEO Series Matrix metadata and CEL files.
options(stringsAsFactors=FALSE, timeout=1800)
set.seed(20261112)

BASE <- "AANAT7223_Ranked_GSEA_Camera_v2"
OUT <- "Supplementary_File_S3_RawCEL_RMA_V25_1"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"RMA_expression"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT,"audit"), recursive=TRUE, showWarnings=FALSE)

required <- c("Biobase","affy","oligo","rat2302cdf","hgu133plus2cdf","pd.ragene.1.0.st.v1")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Required Bioconductor packages missing: ",paste(missing,collapse=", "))

read_series_metadata <- function(gse) {
  f <- file.path(BASE,"rawCEL",paste0(gse,"_matrix"),paste0(gse,"_series_matrix.txt.gz"))
  if(!file.exists(f)) stop("Missing GEO Series Matrix: ",f)
  ln <- readLines(gzfile(f),warn=FALSE)
  get <- function(key) {
    z <- ln[startsWith(ln,paste0(key,"\t"))]
    if(length(z)!=1L) stop("Expected one metadata line for ",key," in ",gse)
    v <- strsplit(z,"\t",fixed=TRUE)[[1]][-1L]
    gsub('"',"",v,fixed=TRUE)
  }
  list(gsm=toupper(get("!Sample_geo_accession")),
       title=get("!Sample_title"),
       platform=get("!Sample_platform_id"))
}

gsm_from_cel <- function(x) sub("^(GSM[0-9]+).*$","\\1",toupper(basename(x)))

dataset_spec <- list(
  GSE23179=list(platform="GPL1355",expected_n=72,selection="All nasal epithelial-cell arrays",method="affy::rma"),
  GSE27263=list(platform="GPL570",expected_n=40,selection="Nasal biopsy arrays only; excludes blood arrays",method="affy::rma"),
  GSE42394=list(platform="GPL6247",expected_n=12,selection="Nose arrays only; excludes white-blood-cell arrays",method="oligo::rma(target='core')")
)

manifest_rows <- list()
summary_rows <- list()
cor_rows <- list()

for(gse in names(dataset_spec)) {
  spec <- dataset_spec[[gse]]
  message("\n===== ",gse," =====")
  meta <- read_series_metadata(gse)
  if(length(meta$gsm)!=length(meta$title)) stop("GSM/title length mismatch: ",gse)
  if(!all(meta$platform==spec$platform)) stop("Platform mismatch for ",gse,": ",paste(unique(meta$platform),collapse=","))
  
  cel_dir <- file.path(BASE,"rawCEL",gse)
  cel <- list.files(cel_dir,pattern="\\.CEL$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
  if(!length(cel)) stop("No decompressed CEL files found for ",gse)
  cel_gsm <- gsm_from_cel(cel)
  if(anyDuplicated(cel_gsm)) stop("Duplicate GSM identifiers in CEL files for ",gse)
  
  if(gse=="GSE23179") {
    selected <- meta$gsm
    reason <- rep("Included: all nasal epithelial-cell samples",length(meta$gsm))
  } else if(gse=="GSE27263") {
    keep <- grepl("^Nasal Biopsy",meta$title,ignore.case=TRUE)
    selected <- meta$gsm[keep]
    reason <- ifelse(keep,"Included: nasal biopsy","Excluded: non-nasal sample")
  } else {
    keep <- grepl("^Nose_",meta$title,ignore.case=TRUE)
    selected <- meta$gsm[keep]
    reason <- ifelse(keep,"Included: nose tissue","Excluded: white blood cells")
  }
  if(length(selected)!=spec$expected_n) stop(gse," selected sample count mismatch: expected ",spec$expected_n,", observed ",length(selected))
  if(!all(selected %in% cel_gsm)) stop(gse," missing selected CEL GSMs: ",paste(setdiff(selected,cel_gsm),collapse=","))
  
  # Preserve full GEO manifest and explicit tissue-selection decisions.
  m <- data.frame(Dataset=gse,GSM=meta$gsm,Title=meta$title,Platform=meta$platform,
                  Selected_for_RMA=meta$gsm %in% selected,
                  Selection_reason=reason,
                  CEL_file=ifelse(meta$gsm %in% cel_gsm,
                                  basename(cel[match(meta$gsm,cel_gsm)]),NA_character_),
                  stringsAsFactors=FALSE)
  manifest_rows[[gse]] <- m
  
  use_cel <- cel[match(selected,cel_gsm)]
  if(anyNA(use_cel) || length(use_cel)!=length(selected)) stop("CEL manifest mapping failed for ",gse)
  
  if(gse %in% c("GSE23179","GSE27263")) {
    suppressPackageStartupMessages(library(affy))
    raw <- affy::ReadAffy(filenames=use_cel)
    Biobase::sampleNames(raw) <- selected
    eset <- affy::rma(raw)
  } else {
    suppressPackageStartupMessages(library(oligo))
    raw <- oligo::read.celfiles(filenames=use_cel)
    Biobase::sampleNames(raw) <- selected
    eset <- oligo::rma(raw,target="core")
  }
  ex <- Biobase::exprs(eset)
  if(!identical(colnames(ex),selected)) {
    if(!setequal(colnames(ex),selected)) stop("RMA output sample names differ from selected GSMs for ",gse)
    ex <- ex[,selected,drop=FALSE]
  }
  if(any(!is.finite(ex))) stop("Non-finite RMA values for ",gse)
  
  cached_file <- file.path(BASE,"RMA_expression",paste0(gse,"_RMA_expression.rds"))
  if(!file.exists(cached_file)) stop("Missing archived comparison RMA matrix: ",cached_file)
  cached <- readRDS(cached_file)
  common_features <- intersect(rownames(ex),rownames(cached))
  common_samples <- intersect(colnames(ex),colnames(cached))
  if(length(common_features)!=nrow(ex) || length(common_features)!=nrow(cached)) {
    message("Feature count differs for ",gse,": new=",nrow(ex),"; cached=",nrow(cached),"; common=",length(common_features))
  }
  if(!setequal(common_samples,selected)) stop("Cached RMA matrix sample IDs differ for ",gse)
  new_aligned <- ex[common_features,selected,drop=FALSE]
  old_aligned <- cached[common_features,selected,drop=FALSE]
  
  for(s in selected) {
    a <- as.numeric(new_aligned[,s]); b <- as.numeric(old_aligned[,s])
    cor_rows[[length(cor_rows)+1L]] <- data.frame(
      Dataset=gse,GSM=s,N_common_features=length(a),
      Pearson_r=cor(a,b,method="pearson"),
      Spearman_rho=cor(a,b,method="spearman"),
      Mean_bias_new_minus_cached=mean(a-b),
      Median_abs_difference=median(abs(a-b)),
      RMSE=sqrt(mean((a-b)^2)),
      stringsAsFactors=FALSE)
  }
  new_rds <- file.path(OUT,"RMA_expression",paste0(gse,"_RMA_expression_reproduced.rds"))
  saveRDS(ex,new_rds,compress="xz")
  summary_rows[[gse]] <- data.frame(
    Dataset=gse,Platform=spec$platform,GEO_samples=length(meta$gsm),
    CEL_files=length(cel),Selected_n=length(selected),
    RMA_features=nrow(ex),RMA_samples=ncol(ex),
    Cached_features=nrow(cached),Cached_samples=ncol(cached),
    Common_features=length(common_features),
    Median_Pearson_r=median(vapply(cor_rows,function(z) if(z$Dataset[1]==gse) z$Pearson_r else NA_real_,numeric(1)),na.rm=TRUE),
    Median_Spearman_rho=median(vapply(cor_rows,function(z) if(z$Dataset[1]==gse) z$Spearman_rho else NA_real_,numeric(1)),na.rm=TRUE),
    stringsAsFactors=FALSE)
  message("RMA completed: ",gse," ",nrow(ex)," x ",ncol(ex),
          "; cached ",nrow(cached)," x ",ncol(cached),
          "; common features ",length(common_features))
}

manifest <- do.call(rbind,manifest_rows)
cor_df <- do.call(rbind,cor_rows)
summary_df <- do.call(rbind,summary_rows)
write.csv(manifest,file.path(OUT,"GEO_CEL_sample_manifest.csv"),row.names=FALSE,na="")
write.csv(cor_df,file.path(OUT,"audit","RMA_vs_cached_per_sample.csv"),row.names=FALSE)
write.csv(summary_df,file.path(OUT,"audit","RMA_reproduction_summary.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
writeLines(c(
  "AANAT7223 raw-CEL/RMA reproducibility — execution record",
  paste("Timestamp UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),
  paste("R:",R.version.string),
  paste("Bioconductor:",as.character(BiocManager::version())),
  "Datasets: GSE23179 (GPL1355, 72 nasal epithelial arrays); GSE27263 (GPL570, 40 nasal-biopsy arrays); GSE42394 (GPL6247, 12 nose arrays).",
  "GSE27263 blood arrays and GSE42394 white-blood-cell arrays are excluded by the pre-specified nasal-tissue selection.",
  "RMA expression matrices are re-derived from the locally archived raw CEL files; no primary limma model is refitted.",
  "Per-sample Pearson/Spearman correlations and expression-difference metrics compare re-derived matrices with the archived RMA matrices used in the pathway analysis."
),file.path(OUT,"README_execution.txt"))
cat("S3_RMA_REPRODUCTION_COMPLETE\n")
print(summary_df)
cat("PER_SAMPLE_QC\n")
print(aggregate(cbind(Pearson_r,Spearman_rho,Mean_bias_new_minus_cached,RMSE)~Dataset,cor_df,
                function(x)c(median=median(x),min=min(x),max=max(x))))
cat("SESSION_INFO_WRITTEN\n")
