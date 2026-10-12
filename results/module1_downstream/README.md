# Module 1 downstream handoff

- `sample_metadata_verified.csv`: all 119 aligned samples, with H mapped to HC as explicitly specified; group totals HC 32, MAP 57, MSAP 20, SAP 10.
- `exploratory_log2_quantile_expression_plus1.tsv.gz`: transformed matrix for exploratory visualization/scoring; source units are unknown, so do not treat as validated quantitative expression for confirmatory pathway scoring.
- `gene_id_symbol_mapping.csv`: source ID/label identity map only; annotation release is unverified.
- `candidate_filtered_source_mean_ge_0.1.tsv.gz` and `candidate_filtered_quantile_normalized.tsv.gz`: reproduce the input filter and normalized values from the uploaded workbook; source units remain unknown.
- Five all-gene exploratory limma tables and all-gene moderated-t ranked lists are available. Ranked lists are not DEG-filtered; use only with the scale caveat above.
- No raw-count matrix was found. Do not use these continuous matrices with DESeq2 or limma-voom.
- No NETosis gene intersection or downstream biological module analysis has been performed.
