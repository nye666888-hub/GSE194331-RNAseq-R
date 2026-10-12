#!/usr/bin/env bash
set -euo pipefail
repo=/workspace/GSE194331-RNAseq-R
cd "$repo"
timestamp=$(date -u +%Y%m%dT%H%M%SZ)
archive="$repo/results/module1_archives/$timestamp"
mkdir -p "$archive/results" "$archive/reports"
for d in module1_audit module1_qc module1_deg module1_severity module1_downstream; do
  if find "results/$d" -maxdepth 1 -type f ! -name original_file_sha256.json -print -quit 2>/dev/null | rg -q .; then
    cp -a "results/$d" "$archive/results/"
  fi
done
if find reports/module1 -maxdepth 1 -type f -print -quit 2>/dev/null | rg -q .; then
  cp -a reports/module1 "$archive/reports/"
fi
mkdir -p /workspace/cloud-setup/logs
log="/workspace/cloud-setup/logs/module1_${timestamp}.log"
exec > >(tee -a "$log") 2>&1
echo "UTC start: $timestamp"
echo 'GitHub source archive: 94abd7389e06991243a8c467c0dcedab5e0ea40f'
python scripts/module1/inventory.py
python scripts/module1/import_original.py
export PATH=/workspace/cloud-setup/r-env/bin:$PATH
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
Rscript --vanilla scripts/module1/run_module1.R
Rscript --vanilla scripts/module1/compare_existing_normalization.R
Rscript --vanilla scripts/module1/run_limma_analysis.R
python scripts/module1/finalize_limma_reports.py
cp "$log" reports/module1/run.log
python scripts/module1/finalize_manifest.py
echo "Log: $log"
git status --short
echo 'MODULE1_EXPLORATORY_LIMMA_COMPLETE; interpret inference as exploratory because source units are unknown.'
