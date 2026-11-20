library(data.table)
base <- "AANAT7223_Ranked_GSEA_Camera_v2/all_feature_statistics"
files <- c(
  GSE42394=file.path(base,"GSE42394_AllFeature_Limma_Conservative95CI.csv.gz"),
  GSE27263=file.path(base,"GSE27263_AllFeature_Limma_Conservative95CI.csv.gz"),
  GSE21477=file.path(base,"GSE21477_AllFeature_Limma_Conservative95CI.csv.gz")
)
for (ds in names(files)) {
  x <- fread(files[[ds]])
  cat("\nDATASET",ds,"ROWS",nrow(x),"\n")
  print(unique(x$Contrast))
  print(x[,.(platform_FDR=sum(BH_FDR<0.05,na.rm=TRUE),
             platform_FDR_1p5=sum(BH_FDR<0.05 & abs(log2FC)>=log2(1.5),na.rm=TRUE)),
          by=Contrast])
}
