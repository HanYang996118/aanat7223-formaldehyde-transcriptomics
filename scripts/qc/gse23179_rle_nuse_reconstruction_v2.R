# Reviewer #2 follow-up: reconstruct per-array RLE/NUSE summaries from the 72 matched raw CEL files
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
for (pkg in c("affy","affyPLM","Biobase","rat2302cdf")) {
  if (!requireNamespace(pkg, quietly=TRUE))
    BiocManager::install(pkg, ask=FALSE, update=FALSE)
}
if (!requireNamespace("data.table", quietly=TRUE))
  install.packages("data.table", repos="https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(affy); library(affyPLM); library(Biobase)
  library(rat2302cdf); library(data.table)
})
base <- "AANAT7223_Ranked_GSEA_Camera_v2"
cel_dir <- file.path(base, "rawCEL", "GSE23179")
matrix_file <- file.path(base, "rawCEL", "GSE23179_matrix", "GSE23179_series_matrix.txt.gz")
stopifnot(dir.exists(cel_dir), file.exists(matrix_file))
cel <- list.files(cel_dir, pattern="\\.CEL$", full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
gsm_from_file <- toupper(regmatches(basename(cel), regexpr("GSM[0-9]+", basename(cel), ignore.case=TRUE)))
cel_map <- data.table(GSM=gsm_from_file, CEL_path=cel)
cel_map <- cel_map[!is.na(GSM) & nzchar(GSM)]
stopifnot(nrow(cel_map)==72L, uniqueN(cel_map$GSM)==72L)
con <- gzfile(matrix_file, "rt")
matrix_lines <- readLines(con, warn=FALSE)
close(con)
acc_line <- matrix_lines[startsWith(matrix_lines,"!Sample_geo_accession")][1]
title_line <- matrix_lines[startsWith(matrix_lines,"!Sample_title")][1]
stopifnot(!is.na(acc_line), !is.na(title_line))
tab <- intToUtf8(9)
acc <- strsplit(acc_line,tab,fixed=TRUE)[[1]][-1]
ttl <- strsplit(title_line,tab,fixed=TRUE)[[1]][-1]
clean_field <- function(x) gsub('"',"",x,fixed=TRUE)
manifest <- data.table(GSM=toupper(clean_field(acc)), title=clean_field(ttl))
stopifnot(nrow(manifest)==72L, uniqueN(manifest$GSM)==72L)
stopifnot(setequal(manifest$GSM, cel_map$GSM))
cel_files <- cel_map$CEL_path[match(manifest$GSM, cel_map$GSM)]
stopifnot(!anyNA(cel_files),all(file.exists(cel_files)))
cat("CEL_MATCH",length(cel_files),"of 72; expected CDF rat2302","\\n")
raw <- affy::ReadAffy(filenames=cel_files,verbose=FALSE)
Biobase::sampleNames(raw) <- manifest$GSM
plm_fit <- affyPLM::fitPLM(raw, background=TRUE, normalize=TRUE, verbosity.level=1)
RLE_matrix <- as.matrix(affyPLM::RLE(plm_fit, type="values"))
NUSE_matrix <- as.matrix(affyPLM::NUSE(plm_fit, type="values"))
stopifnot(nrow(RLE_matrix)==31099L, ncol(RLE_matrix)==72L,
          nrow(NUSE_matrix)==31099L, ncol(NUSE_matrix)==72L)
qc <- data.table(
  GSM=manifest$GSM, title=manifest$title,
  RLE_median=apply(RLE_matrix,2,median,na.rm=TRUE),
  RLE_IQR=apply(RLE_matrix,2,IQR,na.rm=TRUE),
  RLE_MAD=apply(RLE_matrix,2,mad,na.rm=TRUE),
  NUSE_median=apply(NUSE_matrix,2,median,na.rm=TRUE),
  NUSE_IQR=apply(NUSE_matrix,2,IQR,na.rm=TRUE),
  NUSE_MAD=apply(NUSE_matrix,2,mad,na.rm=TRUE)
)
dose_str <- vapply(strsplit(qc$title,"EB_",fixed=TRUE),function(z) strsplit(z[2],"_",fixed=TRUE)[[1]][1],character(1))
duration_str <- vapply(strsplit(qc$title,"EB_",fixed=TRUE),function(z) strsplit(z[2],"_",fixed=TRUE)[[1]][2],character(1))
qc[, dose_ppm := fifelse(dose_str=="07",0.7,as.numeric(dose_str))]
qc[, duration := duration_str]
metrics <- c("RLE_median","RLE_IQR","RLE_MAD","NUSE_median","NUSE_IQR","NUSE_MAD")
limits <- rbindlist(lapply(metrics,function(metric) {
  x <- qc[[metric]]
  med <- median(x,na.rm=TRUE)
  smad <- 1.4826*mad(x,center=med,constant=1,na.rm=TRUE)
  data.table(metric=metric,median=med,scaled_MAD=smad,
    lower_3scaledMAD=med-3*smad,upper_3scaledMAD=med+3*smad)
}))
flag_cols <- character()
for(metric in metrics) {
  flag <- paste0("flag_",metric)
  idx <- match(metric,limits[["metric"]])
  lo <- limits[["lower_3scaledMAD"]][idx]
  hi <- limits[["upper_3scaledMAD"]][idx]
  qc[[flag]] <- qc[[metric]] < lo | qc[[metric]] > hi
  flag_cols <- c(flag_cols,flag)
}
qc[, robust_screen_flag := rowSums(.SD)>0, .SDcols=flag_cols]
qc[, flagged_metrics := apply(.SD,1,function(x) {
  hit <- names(x)[as.logical(x)]
  if(length(hit)) paste(sub("^flag_","",hit),collapse="; ") else ""
}), .SDcols=flag_cols]
group_summary <- qc[,.(n_arrays=.N,
  RLE_median_group_median=median(RLE_median),
  RLE_IQR_group_median=median(RLE_IQR),
  NUSE_median_group_median=median(NUSE_median),
  NUSE_median_min=min(NUSE_median),
  NUSE_median_max=max(NUSE_median),
  NUSE_IQR_group_median=median(NUSE_IQR),
  n_robust_screen_candidates=sum(robust_screen_flag)
),by=.(dose_ppm,duration)][order(dose_ppm,duration)]
out <- file.path(base,"RLE_NUSE_QC")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(qc,file.path(out,"GSE23179_RLE_NUSE_sample_level_review.csv"))
fwrite(qc,file.path(base,"GSE23179_RLE_NUSE_sample_level_review.csv"))
fwrite(limits,file.path(out,"GSE23179_RLE_NUSE_robust_screen_limits.csv"))
fwrite(group_summary,file.path(out,"GSE23179_RLE_NUSE_group_summary.csv"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
saveRDS(RLE_matrix,file.path(out,"GSE23179_RLE_matrix.rds"),compress="xz")
saveRDS(NUSE_matrix,file.path(out,"GSE23179_NUSE_matrix.rds"),compress="xz")
writeLines(c(
 "GSE23179 affyPLM RLE/NUSE reconstruction",
 "72/72 unique GSM-matched CEL files; 31,099 probe sets x 72 arrays.",
 paste("Exploratory flags:",sum(qc$robust_screen_flag)),
 "Flags are manual-review aids only; no array excluded."
),file.path(out,"GSE23179_RLE_NUSE_QC_run_summary.txt"))
zip_path <- file.path(base,"Supplementary_File_S10_GSE23179_RLE_NUSE_QC.zip")
utils::zip(zipfile=zip_path,files=list.files(out,full.names=TRUE),flags="-j")
cat("RLE_NUSE_QC_ZIP",zip_path,"SIZE",file.info(zip_path)$size,"\\n")
cat("RLE_NUSE_COMPLETE arrays=",nrow(qc),
    "flagged=",sum(qc$robust_screen_flag),
    "median_RLE=",median(qc$RLE_median),
    "median_NUSE=",median(qc$NUSE_median),"\\n")
print(limits)
print(group_summary)
print(qc[robust_screen_flag==TRUE,.(GSM,dose_ppm,duration,flagged_metrics)])
