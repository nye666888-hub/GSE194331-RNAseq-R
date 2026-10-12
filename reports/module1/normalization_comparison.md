# Filtering and normalization comparison

The uploaded matrix has 39,376 rows, 119 GSM columns and 3 duplicate-ID rows grouped into 2 duplicated gene identifiers. The existing matrix has 19,645 unique gene rows and the exact same ordered GSM columns. All 19,645 existing IDs occur in the uploaded source.

I used the existing file’s gene set as the comparison filter, then applied average-tie quantile normalization to the unchanged uploaded values. Across every shared gene and all 119 samples, this reproduces the existing matrix: median per-sample Pearson and Spearman correlation 1.000000, maximum per-cell absolute difference 4.511e-10. Candidate output: results/module1_qc/candidate_quantile_normalized_untransformed.tsv.gz. This demonstrates the existing matrix is numerically consistent with quantile-normalizing these source values on its already-filtered gene set. It does not establish the values’ original units, nor validate this as an independent expression filter.

Mean uploaded expression >= 0.1 yields exactly the existing 19,645-gene set (see observed_filter_fingerprint.csv). This is an exact computational fingerprint on these files; it does not confirm the original units or biological rationale. A log2(x+1) followed by quantile normalization did not reproduce the existing matrix. Other filter thresholds are summarized in filter_threshold_sensitivity.csv.

The filtered source and quantile-normalized candidate are saved separately in results/module1_downstream/. The quantile-normalized candidate reproduces the existing matrix within 4.6e-10 maximum absolute cell difference.

The source remains continuous and of unknown provenance. The filtered/quantile-normalized output supports inspection and visual comparison, but a formal DE model still needs the source units and appropriate modeling scale; this output is not raw counts and must not enter DESeq2.
