library(data.table)
dir.create("AANAT7223_Reviewer1_V26_1_Visuals", recursive=TRUE, showWarnings=FALSE)
base <- "AANAT7223_Ranked_GSEA_Camera_v2"
x <- fread(file.path(base,"all_feature_statistics/GSE27263_AllFeature_Limma_Conservative95CI.csv.gz"))
m <- fread(file.path(base,"annotation_mapping/GSE27263_representative_feature_map.csv"))
m <- unique(m[,.(Probe_set=probe_id,gene_symbol)])
z <- merge(x,m,by="Probe_set")
z[,gene_FDR:=p.adjust(P_value,method="BH"),by=Contrast]
out <- z[,.(n_genes=.N,
            gene_FDR=sum(gene_FDR<0.05,na.rm=TRUE),
            gene_FDR_1p5=sum(gene_FDR<0.05 & abs(log2FC)>=log2(1.5),na.rm=TRUE),
            min_gene_FDR=min(gene_FDR,na.rm=TRUE)),by=Contrast]
fk <- z[gene_symbol=="FKBP5",.(Contrast,FKBP5_log2FC=log2FC,FKBP5_P=P_value,FKBP5_gene_FDR=gene_FDR)]
out <- merge(out,fk,by="Contrast",all.x=TRUE)
fwrite(out,"AANAT7223_Reviewer1_V26_1_Visuals/GSE27263_Representative_Gene_FDR_Audit_V26_0.csv")
cat("GSE27263 representative-gene audit from locked all-feature paired limma output\n")
print(out)
