# Remaining plan for confirmatory count-based analysis

- Obtain documented raw integer counts and matching authoritative sample metadata before using DESeq2 or limma-voom. Keep the supplied normalized matrix for concordance checks.
- The exploratory limma contrasts, ordinal score trend, and category global F sensitivity have been run on log2(existing quantile-normalized value + 1). Their source scale is unknown; they are not confirmatory results.
- If valid raw counts become available, fit the predeclared group-only contrasts and severity models, retain all samples, and compare effects only across estimands and scales that are statistically compatible.
- Keep `padj<0.05` and `abs(estimated log2FC)>1` as a descriptive post-test filter unless a thresholded-effect hypothesis is explicitly predeclared. Do not tune cutoffs to reproduce 205 DEGs.
- Verify an authoritative gene annotation release and metadata covariates before adding annotation, age, sex, or batch to future models.
