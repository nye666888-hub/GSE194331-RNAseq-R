# GSE194331 input audit
Normalized candidate has 19645 genes and 119 samples. All sample IDs are GSM-shaped and align one-to-one with local metadata. Compared 19645 shared genes across all 119 exactly aligned samples; both matrices remain continuous with unknown units/provenance.
Group sizes: H/HC 32, MAP 57, MSAP 20, SAP 10. Duplicate gene rows: 0; missing gene IDs: 0; identical expression-vector pairs: 0.
Missing/nonfinite/negative values: 0/0/0. Zero fraction: 0.005919. Integer-like fraction: 0.00596940. Range: 0 to 130629.16.
ID type candidate: gene_symbol_like; no authoritative annotation version
Both continuous expression matrices contain noninteger values and cannot be used as DESeq2 raw counts. They do not identify their expression units or transform.
No preprocessing script, integer raw counts, GEO metadata, RDS object, formal prior DE result or 205-DEG list was found outside explicitly enumerated runtime/cache exclusions. The uploaded source workbook is audited separately.
No values were rounded, renormalized or filtered. Constant rows were omitted only from PCA calculations and remain in the audit/input.
Only supplied-scale descriptive QC was run. Prior cloud-setup/validate.R applied exploratory log2(x+1) for environment testing; it is not evidence of the source scale.
Sample IDs can be checked for duplication; donor-level biological duplication cannot be assessed because donor identifiers are unavailable.
Age/sex/batch variables are absent. Batch-group confounding is unassessable, not disproven. H->HC follows your explicit group definition.
See blocking_report.md and recommended_design.md for formal analysis prerequisites.
