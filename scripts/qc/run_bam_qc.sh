#!/usr/bin/env bash
# Start the BAM QC watcher in the background. Safe to run more than once.
cd "$(dirname "$0")/../.." || exit 1
mkdir -p logs results/qc
if pgrep -f "^python3 scripts/qc/bam_metrics.py --watch" >/dev/null 2>&1; then
  echo "Watcher already running."
else
  nohup python3 scripts/qc/bam_metrics.py --watch >> logs/bam_qc.log 2>&1 &
  echo "Watcher started (PID $!). Log: logs/bam_qc.log"
fi
