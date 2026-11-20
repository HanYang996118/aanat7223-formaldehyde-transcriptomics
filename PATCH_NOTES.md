# Packaging patch notes

## Public-copy-only corrections

The public copies differ from the source-project scripts only in the following static packaging/usability corrections:

1. The helper download_raw_cels() accepts one argument; a call that passed an unsupported second argument was corrected to download_raw_cels(gse).
2. GEOquery and oligo were added to the Bioconductor dependency-installation vector because the script calls both namespaces.
3. The GSE27263 gene-count audit and top-N concordance audit now create their output directory before writing results.

These changes address a function-signature mismatch, missing dependency declarations, output-directory creation, and input-path portability. They do not alter the statistical design, input data, contrasts, gene mapping, pathway definitions, or inferential results. None of the corrected scripts has been represented as freshly executed; run them and verify outputs before describing the archive as fully executable.
