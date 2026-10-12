#!/usr/bin/env python3
import csv, hashlib, json, pathlib
root=pathlib.Path('/workspace/GSE194331-RNAseq-R')
manifest=root/'reports/module1/file_manifest.csv'
rows=[]
for base in [root/'scripts/module1',root/'results/module1_audit',root/'results/module1_qc',
             root/'results/module1_deg',root/'results/module1_severity',root/'results/module1_downstream',root/'reports/module1']:
    for p in sorted(base.rglob('*')):
        if p.is_file() and p!=manifest:
            h=hashlib.sha256()
            with p.open('rb') as f:
                for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
            rows.append({'path':str(p.relative_to(root)),'size_bytes':p.stat().st_size,'sha256':h.hexdigest()})
with manifest.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['path','size_bytes','sha256']);w.writeheader();w.writerows(rows)
baseline=json.loads((root/'results/module1_audit/original_file_sha256.json').read_text())
integrity=[]
for p,expected in baseline.items():
    target=pathlib.Path(p); actual=hashlib.sha256(target.read_bytes()).hexdigest() if target.exists() else None
    integrity.append({'input':str(target.relative_to(root)),'expected_sha256':expected,'actual_sha256':actual,'unchanged':actual==expected})
source=root/'data/module1_inputs/GSE194331.rar'
if source.exists():
    prov=json.loads((root/'results/module1_audit/original_source_provenance.json').read_text())
    actual=hashlib.sha256(source.read_bytes()).hexdigest()
    integrity.append({'input':'data/module1_inputs/GSE194331.rar','expected_sha256':prov['source_sha256'],'actual_sha256':actual,'unchanged':actual==prov['source_sha256']})
with (root/'results/module1_audit/original_integrity_check.csv').open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['input','expected_sha256','actual_sha256','unchanged']);w.writeheader();w.writerows(integrity)
if not all(x['unchanged'] for x in integrity):raise SystemExit('Original-input integrity check failed.')
print(f'Manifested {len(rows)} files; all {len(integrity)} original inputs unchanged.')
