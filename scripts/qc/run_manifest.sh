#!/usr/bin/env bash
cd "$(dirname "$0")/../.." || exit 1
mkdir -p logs results/qc
while IFS= read -r bam; do
  [ -z "$bam" ] && continue
  echo "=== $(date +%H:%M:%S) $bam"
  python3 scripts/qc/bam_metrics.py "$bam"
done < scripts/qc/bam_manifest.txt
echo "=== ALL DONE $(date +%H:%M:%S)"
