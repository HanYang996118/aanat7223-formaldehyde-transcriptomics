# Reproducibility and provenance boundaries

- GSE23179: primary analysis includes 72 arrays. LOAO refits limma after omitting one array while holding the common RMA matrix fixed.
- GSE42394: source-study reproduction applies the documented original signal, fold-change, and ANOVA criteria; it is distinct from the modern limma/BH analyses.
- GSE27263: controlled-human paired nasal-biopsy analyses; the analysed subset is small (8, 7, and 5 paired participants across GEO-recorded lag strata).
- GSE21477: the V26.3.1 manuscript reports recovery/provenance validation of archived outputs (39/39 raw-CEL MD5s; 9/9 historical model summaries), not a new full 39-array rerun. The supplied source-ID-blocked R script is included for transparency only.
- The current script archive is a curated selection from the V26.3/V26.3.1 computational workflow. It is not a claim that every script was freshly executed during repository packaging.
- No raw CEL files, individual-level human data, passwords, tokens, manuscript drafts, or reviewer correspondence are included.
