#!/usr/bin/env python3
"""Inventory analysis inputs without exposing credential files or full matrices."""
import csv, hashlib, json, os, pathlib, zipfile

ROOT = pathlib.Path('/workspace/GSE194331-RNAseq-R')
OUT = ROOT / 'results/module1_audit'
OUT.mkdir(parents=True, exist_ok=True)
rows = []
cache_names = {'r-env', 'mamba', 'cache', 'node_modules', '.git'}
private_names = {'.aws', '.codex', '.agents'}
generated = {'results', 'reports'}
errors = []

def classify(p):
    name = p.name.lower()
    if '归一后' in name:
        return 'expression_candidate', 'Existing expression input', 'unknown_continuous_normalized_candidate'
    if '分组' in name:
        return 'sample_metadata', 'Local sample grouping authority', 'not_applicable'
    if name == 'gse194331.rar':
        return 'uploaded_original_archive', 'Newly uploaded archive; see original_source_provenance.json and original_source_matrix_audit.csv', 'unknown_continuous_source'
    if name == 'gse194331.csv' and p.parent.name == 'extracted':
        return 'ooxml_workbook_mislabeled_csv', 'Actual OOXML workbook imported from GSE194331.rar; see original_source_provenance.json', 'unknown_continuous_source'
    if name == 'gse194331_original_values.tsv.gz':
        return 'derived_original_upload_values', 'Lossless numeric-cell export of uploaded workbook; no scale transformation', 'unknown_continuous_source'
    if name.endswith(('.r', '.rmd', '.py', '.sh')):
        return 'script', 'Inspect provenance; cloud-setup scripts are prior onboarding QC, not original processing', 'not_applicable'
    if name.endswith(('.rds', '.rdata')):
        return 'R_object_candidate', 'Potential counts/metadata/results; requires object inspection', 'unknown'
    if any(t in name for t in ['count', 'tpm', 'fpkm', 'cpm', 'expr', 'deg', 'diff', 'geo', 'gse194331']):
        return 'analysis_candidate', 'Potential analysis input/output; inspect before use', 'unknown'
    if name.endswith(('.csv', '.tsv', '.txt', '.zip', '.gz', '.h5', '.h5ad')):
        return 'tabular_or_archive_candidate', 'Potential data or runtime specification', 'unknown'
    return 'other', 'Documentation, helper, or unrelated file', 'not_applicable'

def add(p, typ=None, purpose=None, scale=None, size=None, sha=''):
    t, u, s = classify(p)
    rows.append(dict(file_path=str(p), type=typ or t, size_bytes=p.stat().st_size if size is None else size,
                     purpose=purpose or u, inferred_expression_scale=scale or s, sha256=sha,
                     scope='workspace_file'))

for base, dirs, files in os.walk('/workspace', onerror=lambda e: errors.append(str(e)), followlinks=False):
    base = pathlib.Path(base)
    for d in list(dirs):
        p = base / d
        if d in private_names or d in cache_names or (base == ROOT and d in generated):
            dirs.remove(d)
            add(p, 'excluded_directory', 'Runtime/cache, protected configuration or current module output; not independent source data', 'not_applicable', 0)
        elif p.is_symlink():
            dirs.remove(d)
            add(p, 'directory_symlink', 'Not followed; inspect manually if relevant', 'unknown', 0)
    for name in sorted(files):
        p = base / name
        try:
            sha = hashlib.sha256(p.read_bytes()).hexdigest() if p.parent == ROOT else ''
            add(p, sha=sha)
            if p.suffix.lower() == '.zip':
                with zipfile.ZipFile(p) as archive:
                    for member in archive.infolist():
                        if member.is_dir():
                            continue
                        virtual = pathlib.Path(str(p) + '::' + member.filename)
                        t, u, s = classify(virtual)
                        rows.append(dict(file_path=str(virtual), type='archive_member_' + t,
                                         size_bytes=member.file_size, purpose=u,
                                         inferred_expression_scale=s, sha256='', scope='zip_member'))
        except Exception as exc:
            errors.append(f'{p}: {exc}')
with (OUT / 'input_inventory.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['file_path','type','size_bytes','purpose','inferred_expression_scale','sha256','scope'])
    writer.writeheader()
    writer.writerows(rows)
(OUT / 'inventory_scope.json').write_text(json.dumps(dict(root='/workspace', errors=errors,
    excluded_runtime_directories=sorted(cache_names), excluded_protected_directories=sorted(private_names),
    notes='Current-module outputs excluded. Prior cloud-setup/results included but are not independent preprocessing evidence.'), ensure_ascii=False, indent=2)+'\n')
if errors:
    raise SystemExit('Inventory errors recorded; inspection is incomplete.')
print(f'Inventory: {len(rows)} entries; no traversal errors.')
