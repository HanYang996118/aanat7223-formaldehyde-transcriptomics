# V26.1 visual regeneration from locked V26.0 source outputs
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

outdir <- "AANAT7223_Reviewer1_V26_1_Visuals"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Figure 5: selected domain terms are derived from the V26.0 complete-family GO output.
go_file <- "metadata/table3_go_display_v26_0.csv"
go <- fread(go_file)
stopifnot(nrow(go) >= 30L,
          all(c("Duration","Concentration_ppm","Direction","Domain","GO_ID","GO_term",
                "Overlap_n","Selected_gene_n","P_value","BH_FDR") %in% names(go)),
          all(go$BH_FDR < 0.05),
          all(go$Concentration_ppm %in% c(6L,10L,15L)))
cat("Figure 5 display rows:",nrow(go),"\n")
print(go[, .N, by=.(Concentration_ppm,Duration,Direction,Domain)][order(Concentration_ppm,Duration,Direction)])

go[, Duration := factor(Duration, levels=c("5 exposure days","4 weeks","13 weeks"))]
go[, Concentration_ppm := factor(Concentration_ppm, levels=c(6,10,15))]
domain_levels <- c("Cell-cycle/proliferative response","DNA damage/repair",
                   "Oxidative/thiol response","Inflammation/immune response",
                   "Ciliary/epithelial structure","Apoptotic response")
go[, Domain := factor(Domain, levels=domain_levels)]
go[, logFDR := -log10(BH_FDR)]
go[, x_id := (as.integer(Concentration_ppm)-1L)*4L +
               match(as.character(Duration), c("5 exposure days","4 weeks","13 weeks"))]
x_breaks <- c(1,2,3,5,6,7,9,10,11)
x_labels <- c("5 d\n6 ppm","4 wk\n6 ppm","13 wk\n6 ppm",
              "5 d\n10 ppm","4 wk\n10 ppm","13 wk\n10 ppm",
              "5 d\n15 ppm","4 wk\n15 ppm","13 wk\n15 ppm")

p5 <- ggplot(go, aes(x=x_id, y=Domain, size=logFDR,
                     colour=Direction, shape=Direction)) +
  geom_vline(xintercept=c(4,8), linetype="22", linewidth=0.45, colour="#AAB6C2") +
  geom_point(alpha=0.92, stroke=0.45, fill="white") +
  scale_x_continuous(breaks=x_breaks, labels=x_labels,
                     expand=expansion(add=c(0.45,0.45))) +
  scale_y_discrete(limits=rev(domain_levels), drop=FALSE) +
  scale_shape_manual(values=c(Up=24,Down=25),
                     labels=c(Up="Upregulated-gene enrichment",Down="Downregulated-gene enrichment")) +
  scale_colour_manual(values=c(Up="#176B87",Down="#C66B36"),
                      labels=c(Up="Upregulated-gene enrichment",Down="Downregulated-gene enrichment")) +
  scale_size_continuous(name=expression(-log[10]("BH-FDR")),
                        range=c(2.6,8.5), breaks=c(2,4,6,8,10),
                        limits=c(0,max(go$logFDR)*1.02)) +
  labs(title="Direction-specific GO-BP over-representation in GSE23179",
       subtitle="One lowest-BH-FDR term per mechanistic domain, exposure condition and direction",
       x="Exposure duration and concentration",
       y=NULL, colour=NULL, shape=NULL,
       caption="All 2,132 eligible GO-BP terms (10–500 mapped background genes) entered BH correction\nfor each non-empty contrast × gene-set list; ≥3 overlaps was a display-only filter.") +
  guides(colour="none", shape=guide_legend(order=1, override.aes=list(size=4.5)),
         size=guide_legend(order=2, override.aes=list(shape=21, colour="#5B6670"))) +
  theme_classic(base_size=12, base_family="Arial") +
  theme(plot.title=element_text(face="bold", size=16, colour="#17324D"),
        plot.subtitle=element_text(size=10.5, colour="#536577"),
        axis.text.x=element_text(size=9.5, colour="#263746", lineheight=0.95),
        axis.text.y=element_text(size=10.5, colour="#263746"),
        axis.title.x=element_text(size=11, face="bold", margin=margin(t=10)),
        panel.grid.major.y=element_line(colour="#E7EDF2", linewidth=0.45),
        panel.grid.major.x=element_blank(),
        legend.position="bottom",
        legend.box="vertical",
        legend.text=element_text(size=9.5),
        plot.caption=element_text(size=8.2, colour="#526170", hjust=0,
                                  margin=margin(t=9)),
        plot.margin=margin(12,18,15,10))
ggsave(file.path(outdir,"Figure5_GO_BP_V26_0.png"), p5,
       width=12.2, height=6.7, dpi=600, bg="white")
ggsave(file.path(outdir,"Figure5_GO_BP_V26_0.pdf"), p5,
       width=12.2, height=6.7, device=cairo_pdf, bg="white")

# Full-data rat scatter for the graphical abstract: all 11,176 matched genes per comparison.
matched_file <- "outputs/data_s15_rat_2ppm_matched_gene_inputs_v26_0.csv.gz"
matched <- fread(matched_file)
stopifnot(nrow(matched)==22352L,
          uniqueN(matched[Comparison=="Early_2ppm_GSE23179_5d_vs_GSE42394_7d",gene_symbol])==11176L,
          uniqueN(matched[Comparison=="Later_2ppm_GSE23179_4w_vs_GSE42394_28d",gene_symbol])==11176L)
matched[,Comparison_label := fifelse(grepl("^Early",Comparison),
                                      "Early comparison: 5 d vs 7 d",
                                      "Later comparison: 4 wk vs 28 d")]
matched[,Comparison_label := factor(Comparison_label,
                                    levels=c("Early comparison: 5 d vs 7 d",
                                             "Later comparison: 4 wk vs 28 d"))]
cor_tab <- matched[, .(
  N=.N,
  rho_moderated_t=cor(Moderated_t_GSE23179,Moderated_t_GSE42394,
                      method="spearman", use="complete.obs"),
  rho_signed_log2FC=cor(log2FC_GSE23179,log2FC_GSE42394,
                        method="spearman", use="complete.obs")
), by=Comparison_label]
fwrite(cor_tab,file.path(outdir,"GraphicalAbstract_Rat_Concordance_Stats_V26_0.csv"))

p_scatter <- ggplot(matched,
                    aes(x=Moderated_t_GSE23179,y=Moderated_t_GSE42394)) +
  geom_point(size=0.42, alpha=0.16, colour="#216B8A") +
  geom_smooth(method="lm", formula=y~x, se=FALSE, linewidth=0.85,
              colour="#C66B36") +
  facet_wrap(~Comparison_label, nrow=1) +
  geom_hline(yintercept=0, linewidth=0.35, colour="#D5DDE4") +
  geom_vline(xintercept=0, linewidth=0.35, colour="#D5DDE4") +
  labs(title="Matched-gene concordance between the two in-vivo rat studies",
       subtitle="All shared one-to-one genes; both contrasts are exposed minus control",
       x="GSE23179 moderated t statistic",
       y="GSE42394 moderated t statistic",
       caption="Spearman correlations are calculated over all 11,176 matched genes in each comparison; signed log2FC sensitivity correlations are reported separately.") +
  theme_classic(base_size=11, base_family="Arial") +
  theme(plot.title=element_text(face="bold",size=15,colour="#17324D"),
        plot.subtitle=element_text(size=10,colour="#536577"),
        strip.text=element_text(face="bold",size=10,colour="#263746"),
        axis.title=element_text(size=10,face="bold"),
        axis.text=element_text(size=9,colour="#263746"),
        plot.caption=element_text(size=8.3,colour="#526170",hjust=0),
        panel.spacing=unit(1.1,"lines"),
        plot.margin=margin(12,12,10,12))
ggsave(file.path(outdir,"GraphicalAbstract_Rat_Scatter_AllMatchedGenes_V26_0.png"),
       p_scatter,width=11.8,height=5.2,dpi=600,bg="white")
ggsave(file.path(outdir,"GraphicalAbstract_Rat_Scatter_AllMatchedGenes_V26_0.pdf"),
       p_scatter,width=11.8,height=5.2,device=cairo_pdf,bg="white")

cat("FIGURE5",file.path(outdir,"Figure5_GO_BP_V26_0.png"),"\n")
cat("FULL_SCATTER",file.path(outdir,"GraphicalAbstract_Rat_Scatter_AllMatchedGenes_V26_0.png"),"\n")
print(cor_tab)
cat("Figure 5 input rows:",nrow(go),"; scatter rows:",nrow(matched),"\n")

options(repr.plot.width=12.2, repr.plot.height=6.7)
print(p5)
