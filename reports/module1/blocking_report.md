# Remaining validation limits

- The provided quantile-normalized matrix has unknown source units. Exploratory limma on log2(QN+1) was run at the user's direction, but these statistics are not validated for confirmatory inference.
- No documented raw integer counts were found; DESeq2 and limma-voom were not run.
- GEO metadata retrieval failed through the environment proxy (HTTP 403). Local group IDs align exactly and counts match the cohort specification supplied by the user; H was mapped to HC from that specification.
- The uploaded workbook's source units and upstream transformation remain unknown despite exact recovery of the existing filter and quantile-normalization output.
- Gene symbols are source labels without a verified external annotation version. Batch, age, sex, and donor duplication are unavailable.
- No original 205-gene table exists locally, preventing a methods-matched overlap comparison.
