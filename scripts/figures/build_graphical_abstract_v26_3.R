# Graphical abstract: fully data-derived V26.0 revision (layout v2)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

outdir <- "AANAT7223_Reviewer1_V26_1_Visuals"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# A. Exact Table 2 gene counts from the V26.0 15-contrast reconciliation.
s2 <- fread("outputs/data_s2_contrast_audit_v26_0.csv")
s2[, Duration := fcase(grepl("^5d_",Contrast),"5 days",
                       grepl("^4w_",Contrast),"4 weeks",
                       grepl("^13w_",Contrast),"13 weeks")]
s2[, Concentration_ppm := as.numeric(sub(".*_([0-9.]+)ppm$","\\1",Contrast))]
s2[, FDR_genes := gene_mapped_FDR_lt_0_05]
s2 <- s2[Concentration_ppm %in% c(0.7,2,6,10,15)]
stopifnot(nrow(s2)==15L)
s2[, Duration := factor(Duration,levels=c("13 weeks","4 weeks","5 days"),labels=c("13w","4w","5d"))]
s2[, Concentration := factor(Concentration_ppm,levels=c(0.7,2,6,10,15),
                             labels=c("0.7","2","6","10","15"))]
pA <- ggplot(s2,aes(x=Concentration,y=Duration,fill=log10(FDR_genes+1))) +
  geom_tile(colour="white",linewidth=1.0,width=0.86,height=0.82) +
  geom_text(aes(label=format(FDR_genes,big.mark=",")),fontface="bold",
            size=2.6,colour="#17324D") +
  scale_fill_gradient(low="#EFF4F7",high="#1F6C8A",
                      name=expression(log[10](n+1))) +
  labs(title="A. Differential-expression burden",
       subtitle="Representative-gene FDR-significant counts",
       x="Formaldehyde (ppm)",y="Exposure duration",
       caption="BH-FDR < 0.05; 14,604-gene background.") +
  theme_minimal(base_size=10,base_family="sans") +
  theme(plot.title=element_text(face="bold",size=12,colour="#17324D"),
        plot.subtitle=element_text(size=9,colour="#526477"),
        axis.title=element_text(face="bold",size=9.5),
        axis.text=element_text(size=8.5,colour="#25384A"),
        panel.grid=element_blank(),
        legend.position="bottom",
        legend.title=element_text(size=9.5),
        legend.text=element_text(size=9.5),
        legend.key.width=unit(0.65,"cm"),
        plot.caption=element_text(size=7.8,colour="#526477",hjust=0),
        plot.margin=margin(7,8,4,7))

# B. Direction-specific GO display from the corrected complete-family V26.0 results.
go <- fread("metadata/table3_go_display_v26_0.csv")
go[,Domain_short:=fcase(
  Domain=="Cell-cycle/proliferative response","Cell cycle",
  Domain=="DNA damage/repair","DNA repair",
  Domain=="Oxidative/thiol response","Oxidative",
  Domain=="Inflammation/immune response","Immune",
  Domain=="Ciliary/epithelial structure","Cilia",
  Domain=="Apoptotic response","Apoptosis",
  default=Domain)]
go[,Domain_short:=factor(Domain_short,
                         levels=c("Cell cycle","DNA repair","Oxidative",
                                  "Immune","Cilia","Apoptosis"))]
go[,Duration:=factor(Duration,levels=c("5 exposure days","4 weeks","13 weeks"),
                     labels=c("5d","4w","13w"))]
go[,Concentration:=factor(Concentration_ppm,levels=c(6,10,15),
                          labels=c("6 ppm","10 ppm","15 ppm"))]
go[,logFDR:=-log10(BH_FDR)]
pB <- ggplot(go,aes(x=Duration,y=Domain_short,size=logFDR,
                    colour=Direction,shape=Direction)) +
  geom_point(alpha=0.92,stroke=0.25) +
  facet_grid(.~Concentration,drop=FALSE) +
  scale_shape_manual(values=c(Up=24,Down=25)) +
  scale_colour_manual(values=c(Up="#176B87",Down="#C66B36")) +
  scale_size_continuous(range=c(1.6,5.0),guide="none") +
  labs(title="B. Direction-specific GO-BP terms",
       subtitle="Complete-family BH results; selected terms",
       x=NULL,y=NULL,colour=NULL,shape=NULL,
       caption="34 selected terms; size = −log10(BH-FDR); 2,132 terms tested; ≥3 overlaps display-only.") +
  theme_minimal(base_size=9,base_family="sans") +
  theme(plot.title=element_text(face="bold",size=12,colour="#17324D"),
        plot.subtitle=element_text(size=9,colour="#526477"),
        strip.text=element_text(face="bold",size=9.5,colour="#25384A"),
        axis.text.x=element_text(size=8,colour="#25384A"),
        axis.text.y=element_text(size=8.5,colour="#25384A"),
        axis.title=element_blank(),
        panel.grid.major.y=element_line(colour="#E5EBF0",linewidth=0.3),
        panel.grid.major.x=element_blank(),
        panel.spacing=unit(0.28,"lines"),
        legend.position="bottom",
        legend.text=element_text(size=9.5),
        legend.key.size=unit(0.25,"cm"),
        plot.caption=element_text(size=7.8,colour="#526477",hjust=0),
        plot.margin=margin(7,7,4,7))

# C. Full-data scatter: 11,176 matched one-to-one genes in each rat comparison.
matched <- fread("outputs/data_s15_rat_2ppm_matched_gene_inputs_v26_0.csv.gz")
stopifnot(nrow(matched)==22352L)
matched[,Comparison_label:=fifelse(grepl("^Early",Comparison),
                                   "Early: 5 d vs 7 d",
                                   "Later: 4 wk vs 28 d")]
matched[,Comparison_label:=factor(Comparison_label,
                                  levels=c("Early: 5 d vs 7 d","Later: 4 wk vs 28 d"))]
cor_tab <- matched[,.(N=.N,
  rho_t=cor(Moderated_t_GSE23179,Moderated_t_GSE42394,method="spearman"),
  rho_fc=cor(log2FC_GSE23179,log2FC_GSE42394,method="spearman")),by=Comparison_label]
cor_tab[,label:=sprintf("Spearman ρ = %.3f\nN = 11,176",rho_t)]
pC <- ggplot(matched,aes(x=Moderated_t_GSE23179,y=Moderated_t_GSE42394)) +
  geom_point(size=0.30,alpha=0.14,colour="#216B8A") +
  geom_smooth(method="lm",formula=y~x,se=FALSE,linewidth=0.75,colour="#C66B36") +
  geom_hline(yintercept=0,linewidth=0.25,colour="#D8E0E7") +
  geom_vline(xintercept=0,linewidth=0.25,colour="#D8E0E7") +
  geom_label(data=cor_tab,aes(x=Inf,y=Inf,label=label),inherit.aes=FALSE,
             hjust=1.04,vjust=1.08,size=2.6,linewidth=0.2,
             fill="white",colour="#25384A") +
  facet_wrap(~Comparison_label,nrow=1,scales="free") +
  labs(title="C. Rat cross-study concordance",
       subtitle="11,176 shared genes per matched 2-ppm comparison",
       x="GSE23179 moderated t",y="GSE42394 moderated t",
       caption=sprintf("Signed log2FC sensitivity: ρ = %.3f (early), %.3f (later).",
                       cor_tab$rho_fc[1],cor_tab$rho_fc[2])) +
  theme_classic(base_size=9,base_family="sans") +
  theme(plot.title=element_text(face="bold",size=12,colour="#17324D"),
        plot.subtitle=element_text(size=9,colour="#526477"),
        strip.text=element_text(face="bold",size=9.5,colour="#25384A"),
        axis.title=element_text(face="bold",size=9.5),
        axis.text=element_text(size=9,colour="#25384A"),
        plot.caption=element_text(size=7.8,colour="#526477",hjust=0),
        panel.spacing=unit(0.35,"lines"),
        plot.margin=margin(7,8,4,7))

# Compact experimental-system cards.
card_grob <- function(title,body,fill="#F5F8FA",accent="#2B6C8A") {
  grobTree(
    roundrectGrob(x=0.5,y=0.5,width=0.98,height=0.96,
                  r=unit(0.08,"in"),gp=gpar(fill=fill,col="#C8D6E0",lwd=0.8)),
    rectGrob(x=0.035,y=0.5,width=0.012,height=0.76,
             gp=gpar(fill=accent,col=NA)),
    textGrob(title,x=0.08,y=0.80,just=c("left","center"),
             gp=gpar(fontface="bold",fontsize=12,col="#17324D")),
    textGrob(body,x=0.08,y=0.45,just=c("left","center"),
             gp=gpar(fontsize=10.5,col="#34495B",lineheight=1.08))
  )
}
cards <- list(
  card_grob("GSE23179  |  IN-VIVO RAT",
            "Whole-body rat inhalation\nControl; 0.7/2/6/10/15 ppm\n5 d / 4 wk / 13 wk"),
  card_grob("GSE42394  |  IN-VIVO RAT",
            "Nose-only rat inhalation\n2 ppm\n7 d / 28 d"),
  card_grob("GSE27263  |  HUMAN",
            "Controlled human inhalation\nBase ≤0.7 ppm; brief peaks ≤0.8 ppm\nPost-program nasal biopsies"),
  card_grob("GSE21477  |  CELLS",
            "Primary human nasal cells\nIn vitro: 20–200 μM\nSingle/repeated schedules"),
  card_grob("OCCUPATIONAL",
            "Anatomy-lab monitoring\nOSHA personal/area records\nTask-related peaks",
            fill="#F7F8F5",accent="#6A8A73")
)
top_cards <- arrangeGrob(grobs=cards,ncol=5,padding=unit(0.06,"in"))

sampling_note <- textGrob(
  "Sampling sites differ: GSE23179, Level II (transitional with some respiratory epithelium; maxilloturbinate unspecified);\nGSE42394, maxilloturbinate (epithelial-subtype proportions unreported).",
  x=0.018,y=0.5,just=c("left","center"),
  gp=gpar(fontsize=9.5,col="#354B5E",lineheight=1.12)
)
sampling_box <- grobTree(roundrectGrob(x=0.5,y=0.5,width=0.995,height=0.95,
                                       r=unit(0.06,"in"),
                                       gp=gpar(fill="#F4F7F9",col="#D6E0E7",lwd=0.7)),
                         sampling_note)

# Human-system findings and interpretive boundaries.
bottom_cards <- list(
  card_grob("GSE27263  |  HUMAN",
            "No probe set or mapped gene met\nBH-FDR < 0.05 in the analysed\npost-program biopsy contrasts.",
            fill="#EEF5F8",accent="#2B6C8A"),
  card_grob("GSE21477  |  HUMAN CELLS",
            "Separate treatment-versus-control\ncontrasts; no direct schedule-effect\ntest.",
            fill="#EEF5F8",accent="#2B6C8A"),
  card_grob("INTERPRETIVE BOUNDARY",
            "Transcriptomic findings do not define\na safe workplace concentration.",
            fill="#F7F8F5",accent="#6A8A73"),
  card_grob("OCCUPATIONAL IMPLICATION",
            "Direct exposure monitoring and\nengineering controls remain essential.",
            fill="#F7F8F5",accent="#6A8A73")
)
bottom <- arrangeGrob(grobs=bottom_cards,ncol=4,padding=unit(0.07,"in"))

title_grob <- textGrob("Concentration-associated but context-sensitive formaldehyde responses",
                       gp=gpar(fontface="bold",fontsize=17,col="#17324D"))
subtitle_grob <- textGrob("Cross-dataset transcriptomic re-analysis relevant to anatomy laboratories",
                          gp=gpar(fontsize=10.5,col="#536B7D"))

rowAB <- arrangeGrob(pA,pB,ncol=2,widths=c(1,1),
                       padding=unit(0.08,"in"))
ga <- arrangeGrob(title_grob,subtitle_grob,top_cards,sampling_box,rowAB,pC,bottom,
                  ncol=1,heights=c(0.38,0.20,1.05,0.38,2.45,1.85,0.85),
                  padding=unit(0.09,"in"))

png(file.path(outdir,"Graphical_Abstract_V26_3_DataDerived.png"),
    width=3900,height=2190,res=300,bg="white")
grid.newpage(); grid.draw(ga); dev.off()
cairo_pdf(file.path(outdir,"Graphical_Abstract_V26_3_DataDerived.pdf"),
          width=13,height=7.3,bg="white")
grid.newpage(); grid.draw(ga); dev.off()

prov <- data.table(
  Panel=c("A","B","C","Human biopsies","Human cells","Rat sampling"),
  Source=c("V26.0 15-contrast Data S2 audit","V26.0 Table3_GO_Display derived from complete 2,132-term GO-BP BH family",
           "V26.0 Data S15: 11,176 matched genes per comparison",
           "GSE27263 source-study analysis","GSE21477 source-study analysis",
           "GSM570620 and GSM1150937 source-level sampling descriptions"),
  Key_values=c("0.7 ppm: 0 genes at all durations; 2 ppm: 5, 7, 1; 6–15 ppm: exact counts shown in heatmap",
               "34 selected term rows; size encodes -log10 BH-FDR; up/down direction preserved",
               sprintf("Moderated-t Spearman rho %.3f and %.3f; signed-log2FC sensitivity %.3f and %.3f",
                       cor_tab$rho_t[1],cor_tab$rho_t[2],cor_tab$rho_fc[1],cor_tab$rho_fc[2]),
               "No platform-wide FDR-significant probe set under analysed conditions; post-program biopsies labelled by GEO lag since 0.7-ppm condition",
               "Primary human nasal epithelial-cell culture; separate treatment-vs-control contrasts",
               "Level II nasal region vs maxilloturbinate; anatomical equivalence not assumed")
)
fwrite(prov,file.path(outdir,"Graphical_Abstract_V26_3_Data_Provenance.csv"))
fwrite(cor_tab,file.path(outdir,"Graphical_Abstract_Rat_Concordance_Stats_V26_0.csv"))
cat("GRAPHICAL_ABSTRACT",file.path(outdir,"Graphical_Abstract_V26_3_DataDerived.png"),"\n")
print(prov)
print(cor_tab)
options(repr.plot.width=13, repr.plot.height=7.3)
grid.newpage(); grid.draw(ga)
