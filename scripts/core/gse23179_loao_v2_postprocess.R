# Post-process QC-linked GSE23179 LOAO v2 outputs (no model refitting)
suppressPackageStartupMessages(library(data.table))
out <- "AANAT7223_GSE23179_LOAO_v2"
res <- file.path(out,"results")
loao <- fread(file.path(res,"GSE23179_LOAO_all_72_arrays_x_15_contrasts.csv"))
qc <- fread("GSE23179_RLE_NUSE_sample_level_review.csv")
qc[,robust_screen_flag:=as.logical(robust_screen_flag)]
stopifnot(nrow(qc)==72L,uniqueN(qc$GSM)==72L,sum(qc$robust_screen_flag)==13L,
          nrow(loao)==1080L,uniqueN(loao$Omitted_GSM)==72L)
role_qc <- loao[,.(N=.N,
 Median_Spearman=median(Spearman_t_vs_primary,na.rm=TRUE),
 Min_Spearman=min(Spearman_t_vs_primary,na.rm=TRUE),
 Median_FDR_Jaccard=if(all(is.na(FDR05_Jaccard))) NA_real_ else median(FDR05_Jaccard,na.rm=TRUE),
 Min_FDR_Jaccard=if(all(is.na(FDR05_Jaccard))) NA_real_ else min(FDR05_Jaccard,na.rm=TRUE),
 Median_P95_abs_delta=median(P95_abs_delta_log2FC,na.rm=TRUE),
 Max_P95_abs_delta=max(P95_abs_delta_log2FC,na.rm=TRUE),
 Max_abs_delta=max(Max_abs_delta_log2FC,na.rm=TRUE),
 Min_Top500_direction=min(Top500_primary_direction_concordance,na.rm=TRUE),
 Max_sign_switch=max(Primary_FDR05_direction_switches,na.rm=TRUE)
),by=.(Omitted_PLM_flag,Omitted_sample_role)][order(Omitted_PLM_flag,Omitted_sample_role)]
fwrite(role_qc,file.path(res,"GSE23179_LOAO_by_PLM_flag_and_role.csv"))
relevant <- loao[Omitted_sample_role %chin% c("exposed","control")]
relevant_qc <- relevant[,.(N=.N,
 Median_Spearman=median(Spearman_t_vs_primary,na.rm=TRUE),
 Min_Spearman=min(Spearman_t_vs_primary,na.rm=TRUE),
 Median_FDR_Jaccard=if(all(is.na(FDR05_Jaccard))) NA_real_ else median(FDR05_Jaccard,na.rm=TRUE),
 Min_FDR_Jaccard=if(all(is.na(FDR05_Jaccard))) NA_real_ else min(FDR05_Jaccard,na.rm=TRUE),
 Median_P95_abs_delta=median(P95_abs_delta_log2FC,na.rm=TRUE),
 Max_P95_abs_delta=max(P95_abs_delta_log2FC,na.rm=TRUE),
 Max_abs_delta=max(Max_abs_delta_log2FC,na.rm=TRUE),
 Min_Top500_direction=min(Top500_primary_direction_concordance,na.rm=TRUE),
 Max_sign_switch=max(Primary_FDR05_direction_switches,na.rm=TRUE)
),by=.(Omitted_PLM_flag,Omitted_sample_role)][order(Omitted_PLM_flag,Omitted_sample_role)]
fwrite(relevant_qc,file.path(res,"GSE23179_LOAO_relevant_omissions_by_PLM_flag.csv"))
fwrite(qc[robust_screen_flag==TRUE][order(GSM)],file.path(res,"GSE23179_13_PLM_flagged_arrays.csv"))
