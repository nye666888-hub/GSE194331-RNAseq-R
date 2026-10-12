Uploaded source: 39376 rows x 119 samples; 3 rows have duplicated identifiers (2 duplicated ID groups).
Existing normalized file: 19645 rows x 119 samples.
Sample columns match exactly: TRUE. Existing target genes found in uploaded source: 19645/19645.
Filter comparison uses the existing matrix feature set only; it reproduces its filtering set and is not an independently chosen expression filter.
Mean expression >= 0.1 retains 19645 unique genes and exactly matches the existing gene set (TRUE); this is an empirical filter fingerprint, not confirmation of units or biological rationale.
Positive-prevalence and other absolute-mean filter profiles are sensitivities because source units are unknown.
Candidate A quantile normalizes the unchanged uploaded values on the common, existing-filter gene set.
Candidate B applies log2(x+1), then quantile normalization, solely as a scale-dependent sensitivity; it is not endorsed absent source units.
Outputs are candidate reprocessing comparisons, not raw counts, a DESeq2 input, or a substitute for provenance.
