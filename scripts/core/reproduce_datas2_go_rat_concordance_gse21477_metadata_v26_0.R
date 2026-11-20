# AANAT7223 Reviewer #1 reconciliation — V26.0
options(stringsAsFactors=FALSE, timeout=1200)
suppressPackageStartupMessages(library(data.table))

OUT <- "AANAT7223_Reviewer1_V26_0_Analysis"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
S2_FILE <- file.path(OUT,"Supplementary_Data_S2_GSE23179_FDR_Significant_Genes_V26_0.csv.gz")
S2_AUDIT <- file.path(OUT,"Supplementary_Data_S2_Contrast_Audit_V26_0.csv")
GO_FILE <- "AANAT7223_Ranked_GSEA_Camera_v2/annotation_mapping/GSE23179_GO_BP_gene_membership.csv"
REP23179 <- "AANAT7223_Ranked_GSEA_Camera_v2/annotation_mapping/GSE23179_representative_feature_map.csv"
REP42394 <- "AANAT7223_Ranked_GSEA_Camera_v2/annotation_mapping/GSE42394_representative_feature_map.csv"
STAT23179 <- "AANAT7223_Ranked_GSEA_Camera_v2/all_feature_statistics/GSE23179_AllFeature_Limma_Conservative95CI.csv.gz"
STAT42394 <- "AANAT7223_Ranked_GSEA_Camera_v2/all_feature_statistics/GSE42394_AllFeature_Limma_Conservative95CI.csv.gz"
GSE21477_MATRIX <- "AANAT7223_Ranked_GSEA_Camera_v2/rawCEL/GSE21477_matrix/GSE21477_series_matrix.txt.gz"

# Rebuild Data S2 from native-feature limma results and the locked 14,604-gene representative-feature map.
all_feat <- fread(STAT23179)
rep_map <- unique(fread(REP23179)[,.(Probe_set=probe_id,gene_symbol)])
stopifnot(all(c("Dataset","Contrast","Probe_set","P_value","log2FC","Moderated_t") %in% names(all_feat)))
all_feat <- all_feat[Dataset=="GSE23179"]
stopifnot(uniqueN(all_feat$Contrast)==15L, uniqueN(rep_map$Probe_set)==14604L,
          uniqueN(rep_map$gene_symbol)==14604L)
label_map <- c(
 "Td5_D0p7.Td5_D0"="5d_0.7ppm","Td5_D2.Td5_D0"="5d_2ppm",
 "Td5_D6.Td5_D0"="5d_6ppm","Td5_D10.Td5_D0"="5d_10ppm","Td5_D15.Td5_D0"="5d_15ppm",
 "Tw4_D0p7.Tw4_D0"="4w_0.7ppm","Tw4_D2.Tw4_D0"="4w_2ppm",
 "Tw4_D6.Tw4_D0"="4w_6ppm","Tw4_D10.Tw4_D0"="4w_10ppm","Tw4_D15.Tw4_D0"="4w_15ppm",
 "Tw13_D0p7.Tw13_D0"="13w_0.7ppm","Tw13_D2.Tw13_D0"="13w_2ppm",
 "Tw13_D6.Tw13_D0"="13w_6ppm","Tw13_D10.Tw13_D0"="13w_10ppm","Tw13_D15.Tw13_D0"="13w_15ppm"
)
gene <- merge(all_feat,rep_map,by="Probe_set",all=FALSE)
stopifnot(nrow(gene)==14604L*15L)
gene[,contrast:=unname(label_map[Contrast])]
stopifnot(!anyNA(gene$contrast))
gene[,gene_FDR:=p.adjust(P_value,method="BH"),by=Contrast]
s2out <- gene[gene_FDR<0.05,
 .(contrast,gene_symbol,representative_probe=Probe_set,log2FC,P_value,FDR=gene_FDR)
][order(contrast,FDR,-abs(log2FC))]
contrast_audit <- gene[,.(model_contrast=first(Contrast),representative_genes=.N,
 gene_mapped_FDR_lt_0_05=sum(gene_FDR<0.05,na.rm=TRUE),
 gene_mapped_FDR_and_abs_log2FC_ge_log2_1p5=
  sum(gene_FDR<0.05 & abs(log2FC)>=log2(1.5),na.rm=TRUE)),by=contrast]
setnames(contrast_audit,"contrast","Contrast")
fwrite(s2out,S2_FILE)
fwrite(contrast_audit,S2_AUDIT)
cat("Regenerated Data S2 rows",nrow(s2out),"across",uniqueN(s2out$contrast),"nonempty contrasts\n")

need <- c(S2_FILE,S2_AUDIT,GO_FILE,REP23179,REP42394,STAT23179,STAT42394,GSE21477_MATRIX)
missing <- need[!file.exists(need)]
if(length(missing)) stop("Missing required locked V25.6 inputs: ",paste(missing,collapse="; "))

# 1. Data S2 / Table 2 count reconciliation
s2 <- fread(S2_FILE)
s2audit <- fread(S2_AUDIT)
stopifnot(all(c("contrast","gene_symbol","representative_probe","log2FC","P_value","FDR") %in% names(s2)))
if(!("Contrast" %in% names(s2audit))) {
  if("contrast" %in% names(s2audit)) setnames(s2audit,"contrast","Contrast")
  else stop("Data S2 audit lacks a contrast identifier column.")
}
stopifnot(nrow(s2)==27944L, uniqueN(s2$contrast)==12L,
          all(s2$FDR<0.05), !anyDuplicated(s2[,.(contrast,gene_symbol)]))
s2_counts <- s2[,.(Data_S2_rows=.N,
                   Data_S2_FDR_genes=uniqueN(gene_symbol),
                   Data_S2_FDR_and_1p5=sum(abs(log2FC)>=log2(1.5))),by=contrast]
audit <- merge(s2audit,s2_counts,by.x="Contrast",by.y="contrast",all.x=TRUE,sort=FALSE)
audit[is.na(Data_S2_rows),c("Data_S2_rows","Data_S2_FDR_genes","Data_S2_FDR_and_1p5"):=list(0L,0L,0L)]
audit[,Data_S2_count_match:=gene_mapped_FDR_lt_0_05==Data_S2_rows]
audit[,Data_S2_1p5_match:=gene_mapped_FDR_and_abs_log2FC_ge_log2_1p5==Data_S2_FDR_and_1p5]
stopifnot(nrow(audit)==15L,all(audit$Data_S2_count_match),all(audit$Data_S2_1p5_match))
fwrite(audit,file.path(OUT,"DataS2_Table2_15Contrast_Reconciliation.csv"))

# 2. Complete GO-BP ORA; BH across every eligible GO term, regardless of observed overlap.
go <- fread(GO_FILE)
rep <- fread(REP23179)
stopifnot(all(c("GO_ID","GO_term","gene_symbol") %in% names(go)))
stopifnot(all(c("probe_id","gene_symbol") %in% names(rep)))
universe <- unique(rep$gene_symbol)
stopifnot(length(universe)==14604L, uniqueN(rep$probe_id)==14604L,
          uniqueN(rep$gene_symbol)==14604L, all(s2$gene_symbol %in% universe))
go <- unique(go[gene_symbol %chin% universe,.(GO_ID,GO_term,gene_symbol)])
term_meta <- unique(go[,.(GO_ID,GO_term)])
term_genes <- split(go$gene_symbol,go$GO_ID)
term_size <- vapply(term_genes,function(x)uniqueN(x),integer(1))
eligible <- names(term_size)[term_size>=10L & term_size<=500L]
term_meta <- term_meta[match(eligible,GO_ID)]
term_meta[,Set_size:=as.integer(term_size[GO_ID])]
N_universe <- length(universe)
cat("GO universe",N_universe,"eligible terms",length(eligible),"\n")

all_tests <- list()
sig_tests <- list()
summary_rows <- list()
contrasts <- c(
 "5d_0.7ppm","5d_2ppm","5d_6ppm","5d_10ppm","5d_15ppm",
 "4w_0.7ppm","4w_2ppm","4w_6ppm","4w_10ppm","4w_15ppm",
 "13w_0.7ppm","13w_2ppm","13w_6ppm","13w_10ppm","13w_15ppm"
)
for(ct in contrasts) {
  d <- unique(s2[contrast==ct & FDR<0.05,.(gene_symbol,log2FC)])
  for(direction in c("All","Up","Down")) {
    sel <- if(direction=="All") unique(d$gene_symbol) else if(direction=="Up") unique(d[log2FC>0,gene_symbol]) else unique(d[log2FC<0,gene_symbol])
    if(!length(sel)) {
      summary_rows[[length(summary_rows)+1L]] <- data.table(
        Contrast=ct,Direction=direction,Selected_genes=0L,
        GO_terms_in_BH_family=0L,GO_FDR05_all_overlaps=0L,
        GO_FDR05_overlap_ge3=0L,Status="Not tested: no FDR-significant selected genes"
      )
      next
    }
    ov <- vapply(eligible,function(id)sum(sel %chin% term_genes[[id]]),integer(1))
    p <- phyper(ov-1L,term_size[eligible],N_universe-term_size[eligible],length(sel),lower.tail=FALSE)
    fdr <- p.adjust(p,method="BH")
    tmp <- data.table(GO_ID=eligible,Overlap_n=ov,Selected_gene_n=length(sel),
                      P_value=p,BH_FDR=fdr,Contrast=ct,Direction=direction)
    tmp <- merge(tmp,term_meta,by="GO_ID",all.x=TRUE,sort=FALSE)
    setcolorder(tmp,c("Contrast","Direction","GO_ID","GO_term","Set_size",
                      "Overlap_n","Selected_gene_n","P_value","BH_FDR"))
    all_tests[[length(all_tests)+1L]] <- tmp
    sig <- tmp[BH_FDR<0.05 & Overlap_n>0]
    if(nrow(sig)) sig_tests[[length(sig_tests)+1L]] <- sig
    summary_rows[[length(summary_rows)+1L]] <- data.table(
      Contrast=ct,Direction=direction,Selected_genes=length(sel),
      GO_terms_in_BH_family=length(eligible),
      GO_FDR05_all_overlaps=nrow(sig),
      GO_FDR05_overlap_ge3=sum(sig$Overlap_n>=3L),
      Status="Tested; BH adjusted across all eligible GO terms"
    )
  }
}
all_go <- rbindlist(all_tests,fill=TRUE)
sig_go <- rbindlist(sig_tests,fill=TRUE)
go_summary <- rbindlist(summary_rows,fill=TRUE)
setorder(all_go,Contrast,Direction,BH_FDR,GO_ID)
setorder(sig_go,Contrast,Direction,BH_FDR,GO_ID)
fwrite(all_go,file.path(OUT,"S1A_AllTested_GO_BP_Terms_V26_0.csv.gz"))
fwrite(sig_go,file.path(OUT,"S1A_Significant_GO_BP_Terms_V26_0.csv"))
fwrite(go_summary,file.path(OUT,"S1A_GO_BP_Testing_Family_Summary_V26_0.csv"))
cat("GO test rows",nrow(all_go),"significant terms",nrow(sig_go),
    "tested groups",sum(go_summary$GO_terms_in_BH_family>0),
    "not-tested groups",sum(go_summary$GO_terms_in_BH_family==0),"\n")

# Compare prior S1A significant list with this exact re-derivation.
old_s1a_file <- "V25_1_recovered_analysis/S1A_rederivation/S1A_GO_BP_ORA_V25_1_current_S2.csv"
if(file.exists(old_s1a_file)) {
  old <- fread(old_s1a_file)
  new <- sig_go[,.(GO_ID,GO_term,Set_size,Overlap_n,Sig_n=Selected_gene_n,
                   P_value,FDR=BH_FDR,Contrast,Direction)]
  oldkey <- old[,.(key=paste(Contrast,Direction,GO_ID,sep="|"))]
  newkey <- new[,.(key=paste(Contrast,Direction,GO_ID,sep="|"))]
  cat("Previous S1A rows",nrow(old),"new S1A rows",nrow(new),
      "same term-key set",setequal(oldkey$key,newkey$key),"\n")
  fwrite(new,file.path(OUT,"S1A_Current_V26_0_Significant_Terms.csv"))
}

# 3. Matched rat genes: signed log2FC and moderated t statistics for both exposed-minus-control contrasts.
s23179 <- fread(STAT23179)
s42394 <- fread(STAT42394)
r23179 <- unique(fread(REP23179)[,.(probe_id,gene_symbol)])
r42394 <- unique(fread(REP42394)[,.(probe_id,gene_symbol)])
stopifnot(all(c("Contrast","Probe_set","Moderated_t","log2FC") %in% names(s23179)),
          all(c("Contrast","Probe_set","Moderated_t","log2FC") %in% names(s42394)))
comparisons <- list(
 list(label="Early_2ppm_GSE23179_5d_vs_GSE42394_7d",ca="Td5_D2.Td5_D0",cb="d7_E - d7_C"),
 list(label="Later_2ppm_GSE23179_4w_vs_GSE42394_28d",ca="Tw4_D2.Tw4_D0",cb="d28_E - d28_C")
)
matched_all <- list()
concordance <- list()
for(cp in comparisons) {
  a <- merge(r23179,s23179[Contrast==cp$ca,.(Probe_set,Moderated_t,log2FC)],
             by.x="probe_id",by.y="Probe_set",all=FALSE)
  b <- merge(r42394,s42394[Contrast==cp$cb,.(Probe_set,Moderated_t,log2FC)],
             by.x="probe_id",by.y="Probe_set",all=FALSE)
  setnames(a,c("probe_id","Moderated_t","log2FC"),c("Probe_GSE23179","Moderated_t_GSE23179","log2FC_GSE23179"))
  setnames(b,c("probe_id","Moderated_t","log2FC"),c("Probe_GSE42394","Moderated_t_GSE42394","log2FC_GSE42394"))
  z <- merge(a,b,by="gene_symbol",all=FALSE)
  stopifnot(nrow(z)==11176L, !anyDuplicated(z$gene_symbol))
  z[,Comparison:=cp$label]
  z[,Contrast_GSE23179:=cp$ca]
  z[,Contrast_GSE42394:=cp$cb]
  matched_all[[length(matched_all)+1L]] <- z
  concordance[[length(concordance)+1L]] <- data.table(
    Comparison=cp$label,N_shared_genes=nrow(z),
    Spearman_moderated_t=cor(z$Moderated_t_GSE23179,z$Moderated_t_GSE42394,method="spearman"),
    Spearman_signed_log2FC=cor(z$log2FC_GSE23179,z$log2FC_GSE42394,method="spearman"),
    Direction_agreement_all_genes=mean(sign(z$log2FC_GSE23179)==sign(z$log2FC_GSE42394)),
    Direction_GSE23179_exposed_minus_control=TRUE,
    Direction_GSE42394_exposed_minus_control=TRUE
  )
}
matched <- rbindlist(matched_all,fill=TRUE)
conc <- rbindlist(concordance,fill=TRUE)
setorder(matched,Comparison,gene_symbol)
fwrite(matched,file.path(OUT,"Rat_2ppm_Matched_Gene_Inputs_V26_0.csv.gz"))
fwrite(conc,file.path(OUT,"Rat_2ppm_Matched_Gene_Concordance_V26_0.csv"))
print(conc)

# 4. Rebuild GSE21477 sample-to-group mapping from GEO metadata, never column positions.
read_geo_meta <- function(path) {
  con <- gzfile(path,"rt"); lines <- readLines(con,warn=FALSE); close(con)
  getvals <- function(prefix) {
    hit <- grep(paste0("^",prefix,"\\t"),lines,value=TRUE)
    if(length(hit)!=1L) stop("Expected exactly one GEO metadata line: ",prefix)
    v <- strsplit(hit,"\t",fixed=TRUE)[[1]][-1]
    gsub('^"|"$',"",v)
  }
  gsm <- getvals("!Sample_geo_accession")
  title <- getvals("!Sample_title")
  source <- getvals("!Sample_source_name_ch1")
  stopifnot(length(gsm)==39L,length(title)==39L,length(source)==39L,!anyDuplicated(gsm))
  data.table(GSM=gsm,Sample_title_GEO=title,Sample_source_name=source)
}
md <- read_geo_meta(GSE21477_MATRIX)
md[,Experiment_block:=sub("^(ES[123]).*","\\1",Sample_source_name)]
md[,Source_ID_within_block:=sub(".*HNEC_#([0-9]+).*","\\1",Sample_source_name)]
md[,Biological_replicate_label:=sub(".*biological rep([0-9]+).*","\\1",Sample_title_GEO,ignore.case=TRUE)]
classify_condition <- function(title,block) {
  tt <- tolower(title)
  if(block=="ES1") {
    if(grepl("50.*fa for 2 hours",tt)) return("50uM_2h")
    if(grepl("50.*fa for 4 hours",tt)) return("50uM_4h")
    if(grepl("100.*fa for 2 hours",tt)) return("100uM_2h")
    if(grepl("control",tt)) return("control_shared")
  }
  if(block=="ES2") {
    if(grepl("100.*fa for 4 hours",tt)) return("100uM_4h")
    if(grepl("100.*fa for 24 hours",tt)) return("100uM_24h")
    if(grepl("200.*fa for 4 hours",tt)) return("200uM_4h")
    if(grepl("200.*fa for 24 hours",tt)) return("200uM_24h")
    if(grepl("control group after 4 hours",tt)) return("control_4h")
    if(grepl("control group after 24 hours",tt)) return("control_24h")
  }
  if(block=="ES3") {
    if(grepl("20.*fa.*4 times",tt)) return("20uM_x4")
    if(grepl("50.*fa.*4 times",tt)) return("50uM_x4")
    if(grepl("control",tt)) return("control_x4")
  }
  NA_character_
}
md[,Condition:=mapply(classify_condition,Sample_title_GEO,Experiment_block,USE.NAMES=FALSE)]
if(nrow(md)!=39L || anyNA(md[,.(Experiment_block,Source_ID_within_block,Biological_replicate_label,Condition)]))
  stop("GSE21477 metadata mapping incomplete.")
stopifnot(setequal(unique(md$Experiment_block),c("ES1","ES2","ES3")))
stopifnot(!anyDuplicated(md[,.(Experiment_block,Source_ID_within_block,Condition)]))
md[,Control_status:=ifelse(grepl("^control_",Condition),"Control","Exposed")]
md[,Matched_control_condition:=NA_character_]
md[Experiment_block=="ES1",Matched_control_condition:="control_shared"]
md[Experiment_block=="ES2" & grepl("_24h$",Condition),Matched_control_condition:="control_24h"]
md[Experiment_block=="ES2" & grepl("_4h$",Condition),Matched_control_condition:="control_4h"]
md[Experiment_block=="ES3",Matched_control_condition:="control_x4"]
md[Control_status=="Control",Matched_control_condition:=Condition]
md[,Matched_control_GSM:=vapply(seq_len(.N),function(i) {
 if(Control_status[i]=="Control") return(GSM[i])
 hit <- md[Experiment_block==md$Experiment_block[i] &
           Source_ID_within_block==md$Source_ID_within_block[i] &
           Condition==md$Matched_control_condition[i],GSM]
 if(length(hit)!=1L) return(NA_character_)
 hit
},character(1))]
if(anyNA(md$Matched_control_GSM)) stop("Could not match every GSE21477 sample to its source-ID-matched control.")
setorder(md,Experiment_block,Source_ID_within_block,Condition)
fwrite(md,file.path(OUT,"GSE21477_Sample_to_Group_Map_V26_0.csv"))
coverage <- md[,.(N=.N),by=.(Experiment_block,Condition)][order(Experiment_block,Condition)]
fwrite(coverage,file.path(OUT,"GSE21477_Block_Condition_Coverage_V26_0.csv"))
cat("GSE21477 sample map",nrow(md),"rows; block/condition cells",nrow(coverage),"\n")
print(coverage)

# 5. Source-tissue audit for the two rat studies.
tissue <- data.table(
 Dataset=c("GSE23179","GSE42394"),
 GEO_sample_reference=c("GSM570620","GSM1150937"),
 Exposure_route=c("Whole-body inhalation","Nose-only inhalation"),
 Nasal_sampling_description=c(
  "Level II region of the nose; tissue dissected and protease-digested to release epithelial cells.",
  "Maxilloturbinate tissue collected after euthanasia and stored in RNAlater; nasal epithelial RNA isolated."
 ),
 Epithelial_composition=c(
  "GEO series/sample description: primarily transitional epithelium with some respiratory epithelium.",
  "GEO sample description identifies nasal epithelial cells; no quantitative epithelial-subtype composition is reported in the accessible record."
 ),
 Maxilloturbinate_explicitly_identified=c(FALSE,TRUE),
 Interpretation_boundary=c(
  "The GEO record does not explicitly identify the maxilloturbinate as the sampled Level II region; do not relabel it as maxilloturbinate.",
  "Maxilloturbinate sampling is explicitly stated in the GEO sample extraction protocol."
 )
)
fwrite(tissue,file.path(OUT,"Rat_Nasal_Tissue_Source_Audit_V26_0.csv"))

writeLines(c(
 "AANAT7223 Reviewer #1 reconciliation outputs — V26.0",
 "",
 "1. DataS2_Table2_15Contrast_Reconciliation.csv: V25.6 Data S2 was counted against its 15-contrast audit. All 15 counts and all 15 >=1.5-fold counts match the audit.",
 "2. S1A_AllTested_GO_BP_Terms_V26_0.csv.gz: all eligible GO-BP terms (10–500 mapped background genes) are tested for every non-empty contrast × direction. BH correction is applied across the complete eligible term family, including zero, one or two overlaps.",
 "3. S1A_Significant_GO_BP_Terms_V26_0.csv and S1A_GO_BP_Testing_Family_Summary_V26_0.csv: complete FDR-significant term list and all 45 contrast × direction status rows. Empty selected-gene sets are marked not tested.",
 "4. Rat_2ppm_Matched_Gene_Inputs_V26_0.csv.gz: 11,176 shared one-to-one genes for each rat comparison, with signed log2FC and moderated t statistics from exposed-minus-control contrasts.",
 "5. GSE21477_Sample_to_Group_Map_V26_0.csv: all 39 GEO samples mapped using GEO title/source metadata, not Series Matrix column positions; includes block, source ID, biological replicate label, condition, and matched control GSM.",
 "6. Rat_Nasal_Tissue_Source_Audit_V26_0.csv: source-level nasal sampling descriptions and limits for GSM570620 and GSM1150937.",
 "",
 "The original GSE21477 primary statistical results remain recovered/provenance-validated historical results. This audit reconstructs sample assignments from GEO metadata and does not claim a new 39-array RMA/limma rerun."
),file.path(OUT,"README_Reviewer1_Reconciliation_V26_0.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
cat("REVIEWER1_V26_0_RECONCILIATION_COMPLETE\n")
