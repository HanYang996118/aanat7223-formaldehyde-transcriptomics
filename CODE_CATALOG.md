# Code catalog and execution guide

This archive contains selected scripts and key result tables from the V26.3 computational workflow. Scripts are grouped by analytical purpose; they are not all standalone because several consume locked intermediate results or journal supplementary archives.

| File | Purpose | Important input / boundary |
|---|---|---|
| scripts/core/reproduce_datas2_go_rat_concordance_gse21477_metadata_v26_0.R | Rebuilds/reconciles Data S2, GO-BP family, rat matched-gene concordance and GSE21477 metadata | Requires locked all-feature limma outputs, representative-feature maps, GO membership and GEO Series Matrix inputs |
| scripts/core/ranked_gsea_camera_v4_intergene_sensitivity.R | Ranked GO-BP GSEA and cameraPR sensitivity across 29 primary contrasts | Requires the all-feature statistics ZIP plus GEO expression matrices/annotation inputs; GSE21477 primary limma statistics are imported, not re-estimated |
| scripts/core/gse23179_loao_v2.R | 72-array × 15-contrast leave-one-array-out limma sensitivity | Requires the common GSE23179 RMA matrix, Series Matrix metadata and locked all-feature archive |
| scripts/core/gse23179_loao_v2_postprocess.R | Summarizes LOAO results by QC flag and sample role | Requires LOAO result CSV and the 72-array RLE/NUSE sample review CSV |
| scripts/preprocessing/rawcel_rma_reproducibility_v25_4.R | Re-derives RMA expression matrices for GSE23179, GSE27263 and GSE42394 | Requires GEO raw CEL files and Series Matrix metadata; does not refit primary limma models |
| scripts/qc/gse23179_rle_nuse_reconstruction_v2.R | Reconstructs probe-level RLE/NUSE metrics for all 72 GSE23179 arrays | Requires 72 matched GSE23179 CEL files and Series Matrix; flags are exploratory, not exclusion rules |
| scripts/source_reproduction/gse42394_source_rule_reproduction_v25_4.R | Reproduces the source-study signal/FC/ANOVA selection rule | Requires deposited GSE42394 Series Matrix; intentionally separate from modern limma/BH inference |
| scripts/secondary/gse21477_sourceid_blocked_limma_v2_transparency_only.R | Contains source-ID-blocked GSE21477 analysis workflow for transparency | The V26.3.1 manuscript reports recovery/provenance validation of archived outputs, not a fresh full 39-array rerun in this revision |
| scripts/figures/build_graphical_abstract_v26_3.R | Rebuilds the V26.3 data-derived graphical abstract | Requires included Data S2 audit and Data S15 inputs, plus the included GO display table |
| scripts/figures/build_visuals.R | Rebuilds Figure 5 and related visual outputs from locked V26.0 tables | Requires the corrected Table 3 GO display output |
| scripts/audit/audit_gse27263_gene_counts.R | Audits representative-gene counts for GSE27263 | Requires the locked GSE27263 all-feature limma table and representative-feature map |
| scripts/audit/audit_table4_allfeature_counts.R | Rechecks platform-wide FDR counts used in Table 4 | Requires all-feature limma result tables for GSE42394, GSE27263 and GSE21477 |
| scripts/audit/audit_table4_topn.R | Recomputes top-N directional concordance from matched-gene statistics | Requires the 11,176-gene matched input table for each rat comparison |

## Recommended execution order

1. Obtain the raw GEO data and the locked supplementary input archives referenced by the scripts.
2. Reconstruct the expression matrices and annotation inputs using the preprocessing workflow where required.
3. Run the Data S2/GO/rat concordance reconciliation against the designated V26.0 inputs.
4. Run the ranked GSEA/cameraPR and LOAO sensitivity workflows only after their expected input files are in place.
5. Run the figure and audit scripts after the locked tables have been regenerated.
6. Compare outputs with the V26.3.1 manuscript and supplementary-file manifest before interpreting any differences.

The outputs/ directory includes the full Data S2 representative-gene table, the complete Data S14 GO-BP test family, the full Data S15 matched-gene input table, and selected RLE/NUSE QC summaries. The supplementary/ directory contains the separate V26.0 Supplementary File S11 ZIP. The four all-feature limma result tables (Supplementary Data S11) are a different supplementary item and are not duplicated in this archive.

Do not infer successful reproduction merely because a script parses or starts. Confirm output dimensions, sample IDs, contrast labels, checksums, and reported summary statistics.
