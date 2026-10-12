# Analysis and method notes

The user requested direct analysis of the supplied normalized expression matrix. We used limma on `log2(QN expression + 1)` as an explicitly exploratory assumption. The underlying unlogged matrix is continuous, unlogged quantile-normalized expression, but its units and processing provenance are not documented. Consequently, the transformation and resulting inferential statistics are not validated for confirmatory differential expression.

No rounded normalized value was represented as a count. DESeq2 and limma-voom were not run because no legal raw-count input was identified. All 19,645 genes were tested; no NETosis-based gene filtering occurred. Binary group contrasts use the verified local cohort labels and full-rank group-only models. The ordinal AP model uses MAP=1, MSAP=2, SAP=3; category-based global F testing is a sensitivity model. SAP n=10 warrants cautious effect-size and confidence-interval interpretation.

BH <0.05 and |estimated log2FC| >1 is only a descriptive post-test DEG filter. Limma tests H0: log2FC=0; no thresholded-effect `treat` test was performed. Moderated-t ranked lists include all finite tested genes and break ties by gene ID. The original 205-DEG result was not available locally for matched-method comparison.
