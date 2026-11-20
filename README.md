# AANAT7223 — Formaldehyde transcriptomics: analysis code

**Code archive version 0.1.0 (release candidate).**  
Associated manuscript: *Cross-dataset transcriptomic re-analysis of formaldehyde reveals concentration-associated but context-sensitive nasal epithelial responses: implications for anatomy laboratories* (Manuscript ID AANAT7223; revision V26.3.1).

## Scope

This repository provides a curated selection of analysis scripts, small provenance/reconciliation tables, and captured software-session records supporting the manuscript revision. It is a code-and-metadata archive; it is **not** a copy of the full raw-data or computational-results archive.

The study re-analyses four public Gene Expression Omnibus (GEO) datasets:

- GSE23179 — rat nasal epithelium, whole-body inhalation
- GSE42394 — rat nasal tissue, nose-only inhalation
- GSE27263 — controlled-human nasal-biopsy study
- GSE21477 — primary human nasal epithelial-cell study

No new human or animal participants were enrolled for this secondary analysis.

## Repository layout

- scripts/core/ — cross-dataset ranked pathway analysis, Data S2/GO/rat concordance reconciliation, and leave-one-array-out (LOAO) analysis
- scripts/preprocessing/ — raw-CEL/RMA preprocessing reproducibility
- scripts/qc/ — GSE23179 RLE/NUSE quality-control reconstruction
- scripts/source_reproduction/ — reproduction of the original GSE42394 source-study criteria
- scripts/secondary/ — GSE21477 source-ID-blocked analysis script, supplied for transparency with the execution boundary described below
- scripts/figures/ — data-derived graphical abstract and visual-generation scripts
- scripts/audit/ — targeted count and table audits
- metadata/ — small, non-participant-level reconciliation and provenance tables
- environment/ — captured R sessionInfo() records from the analysis environments
- outputs/ — key manuscript-facing Data S2, Data S14, Data S15 and RLE/NUSE result tables
- supplementary/ — the actual V26.0 Supplementary File S11 ZIP, included unchanged

See CODE_CATALOG.md for a file-by-file inventory, expected inputs, and execution boundaries. A GitHub Actions workflow performs an R syntax-parse check only; it does not execute the statistical analyses.

## Data access

The key manuscript-facing outputs Data S2, Data S14 and Data S15, together with selected RLE/NUSE QC tables, are included under outputs/. The actual V26.0 Supplementary File S11 archive (ranked GSEA/LOAO/GSE21477 primary provenance bundle) is included unchanged under supplementary/. The distinct Supplementary Data S11 archive containing four platform-wide all-feature limma effect-size/95% CI tables remains a separate item and is not duplicated here. Raw CEL files and GEO Series Matrix records are not redistributed; obtain them from NCBI GEO. Some workflows also require intermediate expression matrices, annotation maps and other supplementary inputs that must be downloaded or regenerated separately.

- GSE23179: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23179
- GSE42394: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42394
- GSE27263: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27263
- GSE21477: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21477

No individual-level participant data, raw CEL files, credentials, API keys, manuscript drafts, or reviewer correspondence are included in this code archive.

## Key analytical rules and interpretation boundaries

1. Native platform-level limma results remain the primary inference; representative-gene and pathway analyses are secondary.
2. The complete GO-BP over-representation family contains 2,132 eligible terms per non-empty contrast–direction list. Benjamini–Hochberg correction includes zero-, one-, and two-overlap terms; the minimum overlap of three is a display-only filter.
3. The corrected Data S14 output contains 72,488 tests across 34 non-empty lists.
4. Rat cross-study comparisons use 11,176 shared one-to-one genes per nominally matched 2-ppm comparison. Weak inverse rank associations are descriptive and are not evidence of biological reversal or independent biological validation.
5. In-vitro micromolar concentrations are not converted to inhaled ppm.
6. GSE23179 LOAO keeps the common 72-array RMA matrix fixed; it is a sample-influence sensitivity analysis, not a re-normalization sensitivity analysis.
7. The 13 RLE/NUSE exploratory flags are manual-review aids; no array was automatically excluded.
8. **GSE21477 execution boundary:** the V26.3.1 revision reports recovery and provenance validation of the archived 39-array source-ID-blocked analysis (39/39 raw-CEL MD5 checks and 9/9 historical model-summary checks). It does not claim a fresh full 39-array RMA/limma rerun for this revision. The corresponding script is included for transparency and must not be cited as evidence that a new full rerun was completed.

## Software environment

The primary cross-dataset workflow was recorded under R 4.4.3 on Debian 12. The GSE21477 provenance/source-ID-blocked script's captured session was R 4.6.2 on Ubuntu 24.04. Exact package versions are recorded in environment/. The workflows use R/Bioconductor packages including limma, fgsea, GEOquery, Biobase, data.table, affy, affyPLM, oligo, and platform-specific annotation packages as required by each script.

The captured session information is descriptive, not a complete lockfile. A fully pinned renv.lock has not been generated in this archive.

## Reproducibility status

This is a curated V26.3/V26.3.1 code-and-metadata selection. It does not claim that every script was executed from a clean environment during packaging. Some scripts require raw GEO files and/or supplementary input archives that are intentionally not redistributed. The manuscript and supplementary materials remain the authoritative record of which analyses were freshly run, recovered, or provenance-validated.

A packaging review corrected one unsupported function argument and added two missing package-installation dependencies in the public ranked GSEA/cameraPR script, redirected figure scripts to the paths inside this archive, and added output-directory creation to two audit scripts. These are static packaging/usability corrections only; they do not represent new statistical analyses or a full rerun. See PATCH_NOTES.md.

## Software authors and affiliations

The software citation metadata follows the manuscript author order:

1. Han Yang — affiliations a, b
2. Rusong Guo — affiliations a, b
3. Yuan Yan — affiliation c
4. Zishuo Lin — affiliation d
5. Mingwang Qiu — affiliations a, b
6. Xiaohui Li — affiliations a, b
7. Kai Zheng — affiliations e, f
8. Zhiyong Fan — affiliations a, b
9. Shiguo Yuan — affiliations e, f
10. Hongwen Liu — affiliations g, h

Affiliations:

- **a** Guangdong Provincial Hospital of Chinese Medicine, Guangzhou 510120, Guangdong, China
- **b** The Second Affiliated Hospital of Guangzhou University of Chinese Medicine, Guangzhou 510120, Guangdong, China
- **c** The Third Affiliated Hospital, Sun Yat-sen University, Guangzhou 510700, China
- **d** School of Traditional Chinese Medicine, Southern Medical University, Guangzhou 510515, Guangdong, China
- **e** Guangdong Provincial Hospital of Chinese Medicine, Hainan Hospital, Guangzhou University of Chinese Medicine, Haikou 570100, Hainan, China
- **f** Hainan Traditional Chinese Medicine Hospital, Hainan Medical University, Haikou 570203, Hainan, China
- **g** Department of Traditional Chinese Orthopedics and Traumatology, Center for Orthopedic Surgery, The Third Affiliated Hospital of Southern Medical University, Guangzhou 510630, Guangdong, China
- **h** Department of Traditional Chinese Medicine, Panzhihua Central Hospital, Panzhihua 617067, Sichuan, China

Han Yang, Rusong Guo, Yuan Yan, and Zishuo Lin contributed equally. Corresponding authors: Zhiyong Fan, Shiguo Yuan, and Hongwen Liu.

## Citation

Please cite the associated manuscript when using these scripts and cite the archived software release once a Zenodo DOI is minted. CITATION.cff provides machine-readable software citation metadata with the ten individual authors in the confirmed manuscript order. Equal contribution: Han Yang, Rusong Guo, Yuan Yan, and Zishuo Lin. Corresponding authors: Zhiyong Fan, Shiguo Yuan, and Hongwen Liu.

## License

The original software code in this repository is released under the MIT License. Third-party GEO-derived data and bundled supplementary materials remain subject to their original source terms. Copyright (c) 2026 Han Yang, Rusong Guo, Yuan Yan, Zishuo Lin, Mingwang Qiu, Xiaohui Li, Kai Zheng, Zhiyong Fan, Shiguo Yuan, and Hongwen Liu. See LICENSE for the complete terms.
