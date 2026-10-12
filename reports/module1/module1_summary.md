# Module 1 summary — exploratory analysis completed with scale limitations

## Inputs, sample audit, and expression scale

- The analyzed expression file is `GSE示例_归一后.zip::GSE示例_归一后.csv`; the supplied group file is `轻中重分组.csv`. Neither original input was modified.
- The expression matrix contains 19,645 genes × 119 samples. The supplied matrix is continuous and numerically identified as unlogged quantile-normalized expression. Its original units and upstream source remain undocumented; it is not raw counts.
- A newly uploaded `GSE194331.rar` contains a mislabeled OOXML workbook with 39,376 gene rows × 119 GSM columns. Mean source expression ≥0.1 recovers the existing unique gene-ID set exactly, and quantile-normalizing those untransformed filtered values reproduces the supplied matrix (maximum absolute difference 4.51e-10). This fingerprints the existing processing but does not establish source units.
- Local labels are H=32, MAP=57, MSAP=20, SAP=10; H was mapped to HC based on your stated cohort definition. All 119 samples were retained. GEO could not be independently queried because the network proxy returned HTTP 403. See the audit outputs for ID, duplicate, missingness, and integrity checks.

## Model and interpretation

- Differential expression used limma 3.68.5 with robust empirical Bayes moderation on `log2(existing quantile-normalized value + 1)`. This transformation is an explicit exploratory modeling assumption because the source units are unknown; it is not asserted to be the original scale.
- Each contrast used a full-rank binary group-only design. No age, sex, or batch covariates were available. The supplied normalized values were neither renormalized nor rounded. All 19,645 genes were tested.
- Reported log2 fold changes are differences on the exploratory transformed scale. P values test H0: log2FC=0; the descriptive DEG rule `BH padj < 0.05 and |estimated log2FC| > 1` is applied after testing, not as an `lfcThreshold` test. These are exploratory results and should not be treated as confirmatory biology until source units are established.
- DESeq2 and limma-voom were not used because no validated raw integer counts were found. The source labels are carried as gene symbols only as unverified identity labels; no external annotation release was confirmed.

## Contrasts

| Comparison | n case vs reference | genes tested | BH < 0.05 | descriptive up | descriptive down |
|---|---:|---:|---:|---:|---:|
| AP vs HC | 87 vs 32 | 19,645 | 10,868 | 112 | 92 |
| MSAP vs MAP | 20 vs 57 | 19,645 | 0 | 0 | 0 |
| SAP vs MAP | 10 vs 57 | 19,645 | 3,266 | 122 | 74 |
| SAP vs MSAP | 10 vs 20 | 19,645 | 0 | 0 | 0 |
| severe vs MAP | 30 vs 57 | 19,645 | 2,927 | 21 | 8 |


The SAP contrasts have only 10 SAP samples. Their estimates and 95% intervals are provided gene by gene, but power and interval width limit interpretation; a nonsignificant result is not evidence of no biological change. Contrast D is exploratory. E is a prespecified severity sensitivity contrast.

## QC and severity

- PCA on the supplied native scale has PC1=68.61% and PC2=15.53%. Because high-magnitude continuous values remain on an unknown source scale, this is descriptive. PCA, distributions, correlation, distance, and group-labeled QC plots are saved in PDF and 300 dpi PNG.
- The screening heuristic flagged 20 candidate samples; these are not established technical outliers. No sample was removed. Age, sex, donor duplication, and technical batches cannot be assessed with the available metadata.
- Among AP samples, the primary ordinal association used MAP=1, MSAP=2, SAP=3 in `~ severity_score`, an 87-sample, rank-2 design (85 residual df); 3,998 genes passed BH<0.05. This tests a constant per-category score slope on the transformed scale and does not prove a biological linear dose response.
- The categorical sensitivity model `~ severity` used a global 2-df moderated F test in a rank-3 design (84 residual df); 2,343 genes passed BH<0.05. All MAP/MSAP/SAP samples were retained.
- Robust empirical Bayes with and without variance-trend moderation was compared for each pairwise contrast. Full per-gene method-stability tables are included.

## Original 205-DEG result and downstream use

- No local result table containing the original 205 DEGs was present, so gene-level overlap and methods-matched comparison cannot be calculated. The descriptive threshold was not tuned to reproduce that count.
- Five complete all-gene limma result tables and moderated-t ranked lists are exported. These can support exploratory pathway review only; validate the source scale and gene annotation before confirmatory GSEA or NETosis claims. The log2(QN+1) matrix and verified four-group metadata are included as exploratory module-two inputs.
- No NETs gene intersection, machine learning, ROC, WGCNA, single-cell, or virtual knockout analysis was performed.

## Reproduction and file index

- Run from the repository root with `bash scripts/module1/run_module1.sh`. The driver archives prior generated outputs, performs inventory/import/audit/QC/normalization fingerprinting, runs limma, refreshes reports, and verifies original input hashes.
- Analysis seed: 194331. R and package versions, parameters, console log, audit files, and SHA-256 file manifest are saved under `reports/module1`, `results/module1_*`, and `scripts/module1`.
- Full generated file inventory: `reports/module1/file_manifest.csv`. Main results are in `results/module1_deg/`; severity results in `results/module1_severity/`; downstream handoff in `results/module1_downstream/`; paired QC plots in `results/module1_qc/`.
- Remaining limitations: undocumented input units/upstream source, unavailable independent GEO annotation retrieval, unverified gene-symbol provenance, and no metadata for batch/age/sex. The exploratory limma results do not resolve those limitations.
