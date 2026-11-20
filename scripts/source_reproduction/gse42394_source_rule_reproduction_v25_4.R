#!/usr/bin/env Rscript
# GSE42394 source-study reproduction, kept separate from the modern limma/BH analysis.
# Operationalization documented in AANAT7223 V25.3 Methods:
# 1) for each duration, use the six nasal arrays (3 control, 3 exposed);
# 2) retain a transcript cluster if mean signal in at least one group is >=
#    the global median of all signal values across the six arrays;
# 3) require absolute fold change >=1.5 (ratio >=1.5 or <=1/1.5);
# 4) require ordinary one-way ANOVA P<0.01 on the deposited linear signal values.
# Expected published counts: 830 (7-day), 42 (28-day).

BASE <- "AANAT7223_Ranked_GSEA_Camera_v2"
OUT <- "Supplementary_File_S3_RawCEL_RMA_V25_1/source_rule"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

f <- file.path(BASE,"rawCEL","GSE42394_matrix","GSE42394_series_matrix.txt.gz")
if(!file.exists(f)) stop("Missing deposited GEO Series Matrix: ",f)
ln <- readLines(gzfile(f),warn=FALSE)

get_meta <- function(key) {
  z <- ln[startsWith(ln,paste0(key,"\t"))]
  if(length(z)!=1L) stop("Expected one metadata row for ",key)
  v <- strsplit(z,"\t",fixed=TRUE)[[1]][-1L]
  gsub('"',"",v,fixed=TRUE)
}
gsm <- get_meta("!Sample_geo_accession")
title <- get_meta("!Sample_title")
st <- grep("^!series_matrix_table_begin",ln)
en <- grep("^!series_matrix_table_end",ln)
if(length(st)!=1L || length(en)!=1L || en<=st) stop("Series Matrix table boundaries not found")
tab <- read.delim(textConnection(ln[(st+1L):(en-1L)]),
                  check.names=FALSE,quote="",comment.char="",
                  stringsAsFactors=FALSE)
if(ncol(tab)!=length(gsm)+1L) stop("Series Matrix expression columns do not match GSM metadata")
expr <- as.matrix(tab[,-1L,drop=FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- as.character(tab[[1L]])
colnames(expr) <- gsm

source_counts <- c("7Day"=830L,"28Day"=42L)
summary_rows <- list()
variant_rows <- list()

for(tp in c("7Day","28Day")) {
  ctrl <- which(grepl(paste0("^Nose_",tp,"Control"),title,ignore.case=TRUE))
  exp <- which(grepl(paste0("^Nose_",tp,"Exposed"),title,ignore.case=TRUE))
  if(length(ctrl)!=3L || length(exp)!=3L) {
    stop(tp," nasal sample selection failed: control=",length(ctrl),"; exposed=",length(exp))
  }
  x <- expr[,c(ctrl,exp),drop=FALSE]
  if(any(!is.finite(x))) stop("Non-finite expression values in ",tp," subset")

  mean_c <- rowMeans(x[,1:3,drop=FALSE])
  mean_e <- rowMeans(x[,4:6,drop=FALSE])
  global_median <- median(x)
  pass_signal <- (mean_c >= global_median) | (mean_e >= global_median)
  fc <- mean_e/mean_c
  log2fc <- log2(fc)

  # Vectorized one-way ANOVA for two groups of three replicates (df1=1, df2=4).
  ssw <- rowSums((x[,1:3,drop=FALSE]-mean_c)^2) +
         rowSums((x[,4:6,drop=FALSE]-mean_e)^2)
  ssb <- 1.5*(mean_e-mean_c)^2
  fstat <- ssb/(ssw/4)
  p_anova <- pf(fstat,df1=1,df2=4,lower.tail=FALSE)

  # Verify the vectorized ANOVA against stats::aov on 100 deterministic rows.
  check_ix <- unique(round(seq(1,nrow(x),length.out=min(100L,nrow(x)))))
  group <- factor(c(rep("Control",3),rep("Exposed",3)),levels=c("Control","Exposed"))
  p_check <- vapply(check_ix,function(i) {
    summary(aov(as.numeric(x[i,])~group))[[1]][["Pr(>F)"]][1]
  },numeric(1))
  if(max(abs(p_anova[check_ix]-p_check),na.rm=TRUE)>1e-8) {
    stop("Vectorized ANOVA failed its aov equivalence check for ",tp)
  }

  responsive <- pass_signal & p_anova<0.01 & (fc>=1.5 | fc<=1/1.5)
  ids <- rownames(x)[which(responsive)]
  published <- unname(source_counts[tp])
  recovered <- length(ids)

  summary_rows[[tp]] <- data.frame(
    Dataset="GSE42394",Timepoint=tp,
    Control_GSM=paste(gsm[ctrl],collapse=";"),
    Exposed_GSM=paste(gsm[exp],collapse=";"),
    Nasal_arrays=6L,Transcript_clusters=nrow(x),
    Global_median_signal=global_median,
    Signal_filter_pass=sum(pass_signal),
    Fold_change_threshold="mean_exposed/mean_control >=1.5 or <=1/1.5",
    ANOVA="ordinary one-way ANOVA on deposited linear signal",
    P_threshold=0.01,Recovered_responsive_clusters=recovered,
    Published_responsive_clusters=published,
    Difference_recovered_minus_published=recovered-published,
    stringsAsFactors=FALSE
  )
  responsive_df <- data.frame(
    ID_REF=ids,Mean_control=mean_c[ids],Mean_exposed=mean_e[ids],
    Fold_change=fc[ids],Log2FC=log2fc[ids],ANOVA_P=p_anova[ids],
    Signal_filter_pass=pass_signal[ids],stringsAsFactors=FALSE
  )
  write.csv(responsive_df,
            file.path(OUT,paste0("GSE42394_",tp,"_Responsive_TranscriptClusters.csv")),
            row.names=FALSE)
  
  # Alternative scale choices are retained as a transparent sensitivity audit.
  xlog <- log2(pmax(x,1e-8))
  mc_log <- rowMeans(xlog[,1:3,drop=FALSE])
  me_log <- rowMeans(xlog[,4:6,drop=FALSE])
  ssw_log <- rowSums((xlog[,1:3,drop=FALSE]-mc_log)^2) +
              rowSums((xlog[,4:6,drop=FALSE]-me_log)^2)
  ssb_log <- 1.5*(me_log-mc_log)^2
  p_log <- pf(ssb_log/(ssw_log/4),df1=1,df2=4,lower.tail=FALSE)
  variant_rows[[tp]] <- data.frame(
    Timepoint=tp,Raw_signal_ANOVA_FC_ratio_1p5=recovered,
    Log2_signal_ANOVA_FC_ratio_1p5=sum(pass_signal&p_log<0.01&(fc>=1.5|fc<=1/1.5),na.rm=TRUE),
    Raw_signal_ANOVA_abs_log2FC_1p5=sum(pass_signal&p_anova<0.01&abs(log2fc)>=1.5,na.rm=TRUE),
    Log2_signal_ANOVA_abs_log2FC_1p5=sum(pass_signal&p_log<0.01&abs(log2fc)>=1.5,na.rm=TRUE),
    stringsAsFactors=FALSE
  )
}

summary_df <- do.call(rbind,summary_rows)
variant_df <- do.call(rbind,variant_rows)
write.csv(summary_df,file.path(OUT,"GSE42394_SourceRule_Reproduction_Summary.csv"),row.names=FALSE)
write.csv(variant_df,file.path(OUT,"GSE42394_SourceRule_Variant_Audit.csv"),row.names=FALSE)
if(!identical(as.integer(summary_df$Recovered_responsive_clusters),c(829L,42L))) {
  stop("Source-rule counts differ from the reconciled expected reproduction (829, 42).")
}
writeLines(c(
  "GSE42394 source-study reproduction — V25.4",
  "Data source: deposited GSE42394 Series Matrix; this is intentionally separate from raw-CEL/RMA + limma/BH primary analysis.",
  "Tissue: nasal respiratory mucosa only; 7-day and 28-day duration-specific comparisons, each with 3 control and 3 exposed arrays.",
  "Low-signal filter: retain a transcript cluster if the mean signal in at least one exposure group is not below the global median of all signal values across the six nasal arrays.",
  "Differential-expression rule: exposed/control mean signal ratio >=1.5 or <=1/1.5, plus ordinary one-way ANOVA P<0.01 on deposited linear signal values.",
  "The vectorized F test was checked against stats::aov for 100 deterministic transcript clusters in each timepoint.",
  paste("Recovered counts:",paste(summary_df$Timepoint,summary_df$Recovered_responsive_clusters,collapse="; ")),
  "Published reference counts: 830 (7-day), 42 (28-day).",
  "Interpretation: 829/830 transcript clusters recovered at 7 days and 42/42 at 28 days; the 7-day difference is one cluster and is not silently forced to match.",
  "Source: Rager et al., 2014, Formaldehyde-Associated Changes in microRNAs: Tissue and Temporal Specificity in the Rat Nose, White Blood Cells, and Bone Marrow."
),file.path(OUT,"README_source_rule.txt"))
cat("SOURCE_RULE_REPRODUCTION_COMPLETE\n")
print(summary_df)
print(variant_df)
