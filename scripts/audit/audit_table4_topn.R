library(data.table)
dir.create("AANAT7223_Reviewer1_V26_1_Visuals", recursive=TRUE, showWarnings=FALSE)
x <- fread("AANAT7223_Reviewer1_V26_0_Analysis/Rat_2ppm_Matched_Gene_Inputs_V26_0.csv.gz")
stopifnot(nrow(x)==22352L, uniqueN(x$gene_symbol)==11176L)
out <- list()
for (cmp in unique(x$Comparison)) {
  d <- x[Comparison==cmp]
  stopifnot(nrow(d)==11176L, uniqueN(d$gene_symbol)==11176L)
  a <- order(-abs(d$Moderated_t_GSE23179), d$gene_symbol)
  b <- order(-abs(d$Moderated_t_GSE42394), d$gene_symbol)
  for (n in c(100L,250L,500L,1000L)) {
    ia <- a[seq_len(n)]; ib <- b[seq_len(n)]
    dir_ab <- mean(sign(d$Moderated_t_GSE23179[ia]) == sign(d$Moderated_t_GSE42394[ia]))
    dir_ba <- mean(sign(d$Moderated_t_GSE42394[ib]) == sign(d$Moderated_t_GSE23179[ib]))
    ov <- length(intersect(d$gene_symbol[ia],d$gene_symbol[ib]))/n
    out[[length(out)+1L]] <- data.table(
      Comparison=cmp,N_top=n,Direction_A_to_B=dir_ab,
      Direction_B_to_A=dir_ba,TopN_overlap=ov,
      Ranking_statistic="absolute moderated t",
      Tie_break="gene_symbol ascending"
    )
  }
}
res <- rbindlist(out)
fwrite(res,"AANAT7223_Reviewer1_V26_1_Visuals/Rat_2ppm_TopN_Direction_Concordance_V26_0.csv")
cat("TOP_N_RECOMPUTED_FROM_DATA_S15\n")
print(res)
cat("RANGES\n")
print(res[,.(A_to_B_min=min(Direction_A_to_B),A_to_B_max=max(Direction_A_to_B),
             B_to_A_min=min(Direction_B_to_A),B_to_A_max=max(Direction_B_to_A)),
          by=Comparison])
