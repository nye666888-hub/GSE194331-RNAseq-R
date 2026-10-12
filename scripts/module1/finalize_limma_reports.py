#!/usr/bin/env python3
"""Refresh module 1 status and reports from the completed exploratory limma run."""
import csv
import pathlib

ROOT = pathlib.Path('/workspace/GSE194331-RNAseq-R')


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def n(x):
    return f'{int(x):,}'


deg = read_csv(ROOT / 'results/module1_deg/limma_DEG_summary.csv')
severity = read_csv(ROOT / 'results/module1_severity/severity_model_summary.csv')
assert len(deg) == 5 and len(severity) == 2
by_name = {row['comparison']: row for row in deg}
ord_sig = next(int(x['significant_BH_lt_0_05']) for x in severity if x['model'] == 'ordinal_score')
cat_sig = next(int(x['significant_BH_lt_0_05']) for x in severity if x['model'] == 'AP_category_global_F')

status_rows = [
    ('inventory', 'COMPLETE'),
    ('sample_id_alignment', 'PASS'),
    ('expression_scale_confirmation', 'SOURCE_UNITS_UNKNOWN_EXPLORATORY_MODEL_ONLY'),
    ('GEO_group_confirmation', 'BLOCKED_NETWORK; local group metadata cross-checked to user-specified counts'),
    ('candidate_design_rank', 'PASS_ALGEBRA_ONLY; binary contrasts and severity designs full rank'),
    ('descriptive_QC', 'COMPLETE_WITH_SCALE_LIMITATIONS'),
    ('exploratory_limma_DE_A_to_E', 'COMPLETE_LOG2_QN_PLUS1_ASSUMPTION'),
    ('exploratory_AP_trend', 'COMPLETE_ORDINAL_AND_CATEGORICAL_SENSITIVITY'),
    ('method_stability', 'COMPLETE_EBAYES_TREND_SENSITIVITY'),
    ('downstream_expression', 'EXPLORATORY_LOG2_QN_PLUS1_MATRIX_EXPORTED'),
    ('GSEA_rank_export', 'EXPLORATORY_MODERATED_T_RANKS_EXPORTED; validate scale before confirmatory use'),
    ('DESeq2_and_voom', 'NOT_RUN_NO_VALIDATED_RAW_COUNTS'),
]
write_csv(ROOT / 'reports/module1/run_status.csv',
          [{'step': a, 'status': b} for a, b in status_rows], ['step', 'status'])

comparison_rows = []
for row in deg:
    comparison_rows.append({
        'comparison': row['comparison'], 'contrast': row['contrast'],
        'n_case': row['n_case'], 'n_reference': row['n_reference'],
        'genes_tested': row['genes_tested'],
        'padj_lt_0_05': row['genes_padj_lt_0_05'],
        'descriptive_DEG_up': row['descriptive_DEG_up'],
        'descriptive_DEG_down': row['descriptive_DEG_down'],
        'status': 'EXPLORATORY_LIMMA_LOG2_QN_PLUS1',
    })
write_csv(ROOT / 'results/module1_deg/comparison_status.csv', comparison_rows)

refs_path = ROOT / 'results/module1_downstream/input_references.csv'
refs = read_csv(refs_path)
ref_status = {
    'supplied_expression_matrix': 'available; unlogged quantile-normalized; source units unknown; used only after explicit exploratory log2(value+1) assumption',
    'raw_counts': 'NOT_FOUND; DESeq2/voom not run',
    'validated_transformed_matrix': 'NOT_VALIDATED; exploratory log2(QN+1) matrix exported with scale caveat',
    'GSEA_ranked_lists': 'EXPLORATORY_MODERATED_T_RANKS_EXPORTED; all tested genes; scale assumption applies',
    'DE_results': 'EXPLORATORY_LIMMA_COMPLETE; not confirmatory because source units unknown',
}
for row in refs:
    if row['item'] in ref_status:
        row['status'] = ref_status[row['item']]
write_csv(refs_path, refs)

summary = f'''# Module 1 summary — exploratory analysis completed with scale limitations

## Inputs, sample audit, and expression scale

- The analyzed expression file is `GSE示例_归一后.zip::GSE示例_归一后.csv`; the supplied group file is `轻中重分组.csv`. Neither original input was modified.
- The expression matrix contains 19,645 genes × 119 samples. The supplied matrix is continuous and numerically identified as unlogged quantile-normalized expression. Its original units and upstream source remain undocumented; it is not raw counts.
- A newly uploaded `GSE194331.rar` contains a mislabeled OOXML workbook with 39,376 gene rows × 119 GSM columns. Mean source expression ≥0.1 recovers the existing unique gene-ID set exactly, and quantile-normalizing those untransformed filtered values reproduces the supplied matrix (maximum absolute difference 4.51e-10). This fingerprints the existing processing but does not establish source units.
- Local labels are H=32, MAP=57, MSAP=20, SAP=10; H was mapped to HC based on your stated cohort definition. All 119 samples were retained. GEO could not be independently queried because the network proxy returned HTTP 403. See the audit outputs for ID, duplicate, missingness, and integrity checks.

## Model and interpretation

- Differential expression used limma {read_csv(ROOT / 'results/module1_deg/limma_package_versions.csv')[1]['version']} with robust empirical Bayes moderation on `log2(existing quantile-normalized value + 1)`. This transformation is an explicit exploratory modeling assumption because the source units are unknown; it is not asserted to be the original scale.
- Each contrast used a full-rank binary group-only design. No age, sex, or batch covariates were available. The supplied normalized values were neither renormalized nor rounded. All 19,645 genes were tested.
- Reported log2 fold changes are differences on the exploratory transformed scale. P values test H0: log2FC=0; the descriptive DEG rule `BH padj < 0.05 and |estimated log2FC| > 1` is applied after testing, not as an `lfcThreshold` test. These are exploratory results and should not be treated as confirmatory biology until source units are established.
- DESeq2 and limma-voom were not used because no validated raw integer counts were found. The source labels are carried as gene symbols only as unverified identity labels; no external annotation release was confirmed.

## Contrasts

| Comparison | n case vs reference | genes tested | BH < 0.05 | descriptive up | descriptive down |
|---|---:|---:|---:|---:|---:|
'''
for key in ['A_AP_vs_HC', 'B_MSAP_vs_MAP', 'C_SAP_vs_MAP', 'D_SAP_vs_MSAP', 'E_MSAP_SAP_vs_MAP']:
    r = by_name[key]
    summary += (f"| {r['contrast']} | {r['n_case']} vs {r['n_reference']} | {n(r['genes_tested'])} | "
                f"{n(r['genes_padj_lt_0_05'])} | {n(r['descriptive_DEG_up'])} | {n(r['descriptive_DEG_down'])} |\n")
summary += f'''

The SAP contrasts have only 10 SAP samples. Their estimates and 95% intervals are provided gene by gene, but power and interval width limit interpretation; a nonsignificant result is not evidence of no biological change. Contrast D is exploratory. E is a prespecified severity sensitivity contrast.

## QC and severity

- PCA on the supplied native scale has PC1=68.61% and PC2=15.53%. Because high-magnitude continuous values remain on an unknown source scale, this is descriptive. PCA, distributions, correlation, distance, and group-labeled QC plots are saved in PDF and 300 dpi PNG.
- The screening heuristic flagged 20 candidate samples; these are not established technical outliers. No sample was removed. Age, sex, donor duplication, and technical batches cannot be assessed with the available metadata.
- Among AP samples, the primary ordinal association used MAP=1, MSAP=2, SAP=3 in `~ severity_score`, an 87-sample, rank-2 design (85 residual df); {n(ord_sig)} genes passed BH<0.05. This tests a constant per-category score slope on the transformed scale and does not prove a biological linear dose response.
- The categorical sensitivity model `~ severity` used a global 2-df moderated F test in a rank-3 design (84 residual df); {n(cat_sig)} genes passed BH<0.05. All MAP/MSAP/SAP samples were retained.
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
'''
(ROOT / 'reports/module1/module1_summary.md').write_text(summary, encoding='utf-8')

notes = '''# Analysis and method notes

The user requested direct analysis of the supplied normalized expression matrix. We used limma on `log2(QN expression + 1)` as an explicitly exploratory assumption. The underlying unlogged matrix is continuous, unlogged quantile-normalized expression, but its units and processing provenance are not documented. Consequently, the transformation and resulting inferential statistics are not validated for confirmatory differential expression.

No rounded normalized value was represented as a count. DESeq2 and limma-voom were not run because no legal raw-count input was identified. All 19,645 genes were tested; no NETosis-based gene filtering occurred. Binary group contrasts use the verified local cohort labels and full-rank group-only models. The ordinal AP model uses MAP=1, MSAP=2, SAP=3; category-based global F testing is a sensitivity model. SAP n=10 warrants cautious effect-size and confidence-interval interpretation.

BH <0.05 and |estimated log2FC| >1 is only a descriptive post-test DEG filter. Limma tests H0: log2FC=0; no thresholded-effect `treat` test was performed. Moderated-t ranked lists include all finite tested genes and break ties by gene ID. The original 205-DEG result was not available locally for matched-method comparison.
'''
(ROOT / 'reports/module1/limma_analysis_notes.md').write_text(notes, encoding='utf-8')

blocked = '''# Remaining validation limits

- The provided quantile-normalized matrix has unknown source units. Exploratory limma on log2(QN+1) was run at the user's direction, but these statistics are not validated for confirmatory inference.
- No documented raw integer counts were found; DESeq2 and limma-voom were not run.
- GEO metadata retrieval failed through the environment proxy (HTTP 403). Local group IDs align exactly and counts match the cohort specification supplied by the user; H was mapped to HC from that specification.
- The uploaded workbook's source units and upstream transformation remain unknown despite exact recovery of the existing filter and quantile-normalization output.
- Gene symbols are source labels without a verified external annotation version. Batch, age, sex, and donor duplication are unavailable.
- No original 205-gene table exists locally, preventing a methods-matched overlap comparison.
'''
(ROOT / 'reports/module1/blocking_report.md').write_text(blocked, encoding='utf-8')

count_limit = '''# Count-based inference not performed

This file concerns count-based methods only. No documented raw integer-count input was identified, so DESeq2 and limma-voom were not run. At the user's direction, exploratory limma results are available under `results/module1_deg/`, using log2(existing quantile-normalized expression + 1). The source units remain unknown, so that analysis is explicitly exploratory and not a substitute for count-based inference.
'''
(ROOT / 'results/module1_deg/BLOCKED.md').write_text(count_limit, encoding='utf-8')
(ROOT / 'results/module1_severity/BLOCKED.md').write_text(count_limit, encoding='utf-8')

downstream = f'''# Module 1 downstream handoff

- `sample_metadata_verified.csv`: all 119 aligned samples, with H mapped to HC as explicitly specified; group totals HC 32, MAP 57, MSAP 20, SAP 10.
- `exploratory_log2_quantile_expression_plus1.tsv.gz`: transformed matrix for exploratory visualization/scoring; source units are unknown, so do not treat as validated quantitative expression for confirmatory pathway scoring.
- `gene_id_symbol_mapping.csv`: source ID/label identity map only; annotation release is unverified.
- `candidate_filtered_source_mean_ge_0.1.tsv.gz` and `candidate_filtered_quantile_normalized.tsv.gz`: reproduce the input filter and normalized values from the uploaded workbook; source units remain unknown.
- Five all-gene exploratory limma tables and all-gene moderated-t ranked lists are available. Ranked lists are not DEG-filtered; use only with the scale caveat above.
- No raw-count matrix was found. Do not use these continuous matrices with DESeq2 or limma-voom.
- No NETosis gene intersection or downstream biological module analysis has been performed.
'''
(ROOT / 'results/module1_downstream/README.md').write_text(downstream, encoding='utf-8')

plan = '''# Remaining plan for confirmatory count-based analysis

- Obtain documented raw integer counts and matching authoritative sample metadata before using DESeq2 or limma-voom. Keep the supplied normalized matrix for concordance checks.
- The exploratory limma contrasts, ordinal score trend, and category global F sensitivity have been run on log2(existing quantile-normalized value + 1). Their source scale is unknown; they are not confirmatory results.
- If valid raw counts become available, fit the predeclared group-only contrasts and severity models, retain all samples, and compare effects only across estimands and scales that are statistically compatible.
- Keep `padj<0.05` and `abs(estimated log2FC)>1` as a descriptive post-test filter unless a thresholded-effect hypothesis is explicitly predeclared. Do not tune cutoffs to reproduce 205 DEGs.
- Verify an authoritative gene annotation release and metadata covariates before adding annotation, age, sex, or batch to future models.
'''
(ROOT / 'reports/module1/analysis_plan_not_executed.md').write_text(plan, encoding='utf-8')

print('Updated module 1 status and reports from generated limma and severity tables.')
