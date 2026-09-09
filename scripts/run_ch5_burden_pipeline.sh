#!/usr/bin/env bash
# Regenerate every Chapter 5 burden table and figure on the primary-contig
# set. Run from the repository root with the conda environment active.
#
#   bash scripts/run_ch5_burden_pipeline.sh
#
# Steps
#   0  preserve the pre-restriction CSV as *_allcontigs.csv (once)
#   1  burden CSV (primary contigs)          scripts/ch5_variant_burden.py
#   2  burden Word tables                    scripts/make_ch5_variant_burden_table.py
#   3  PASS burden figure (single panel)     scripts/plot_ch5_pass_burden.py
#   4  filter-class proportions figure       scripts/plot_ch5_filter_classes.py
#   5  SIFT Fisher tests + table + forest    scripts/ch5_sift_fisher_test.py, make_ch5_sift_fisher_table.py, plot_ch5_sift_fisher_forest.py
#   6  filter cascade, all VCFs under data/  scripts/qc/filter_cascade.py --all --force
#   7  cascade Word tables (MCF-7)           scripts/make_ch5_filter_cascade_table.py
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== 0  preserve all-contig CSV"
if [ -f outputs/ch5_variant_burden.csv ] && [ ! -f outputs/ch5_variant_burden_allcontigs.csv ]; then
  if ! grep -q "contig_set" outputs/ch5_variant_burden.csv; then
    cp outputs/ch5_variant_burden.csv outputs/ch5_variant_burden_allcontigs.csv
    echo "   copied pre-restriction CSV to outputs/ch5_variant_burden_allcontigs.csv"
  fi
fi

echo "== 1  burden CSV (primary contigs)"
python3 scripts/ch5_variant_burden.py

echo "== 2  burden Word tables"
python3 scripts/make_ch5_variant_burden_table.py

echo "== 3  PASS burden figure"
python3 scripts/plot_ch5_pass_burden.py

echo "== 4  filter-class proportions figure"
python3 scripts/plot_ch5_filter_classes.py

echo "== 5  SIFT Fisher tests"
python3 scripts/ch5_sift_fisher_test.py
[ -f scripts/make_ch5_sift_fisher_table.py ] && python3 scripts/make_ch5_sift_fisher_table.py
[ -f scripts/plot_ch5_sift_fisher_forest.py ] && python3 scripts/plot_ch5_sift_fisher_forest.py

echo "== 6  filter cascade (primary contigs, all VCFs under data/)"
python3 scripts/qc/filter_cascade.py --all --force
python3 scripts/qc/filter_cascade.py --table

echo "== 7  cascade Word tables"
python3 scripts/make_ch5_filter_cascade_table.py

echo
echo "Done. Deliverables:"
echo "  outputs/ch5_variant_burden.csv"
echo "  outputs/Table_5_X_variant_burden_{full,condensed}.docx"
echo "  outputs/Table_5_X_filter_cascade.docx"
echo "  figures/mcf7/pass_burden_mcf7.{png,pdf}"
echo "  figures/mcf7/filter_class_proportions_mcf7.{png,pdf}"
echo "  figures/mcf7/sift_fisher_forest_mcf7.{png,pdf}"
echo "  results/cascade_primary/cascade_summary.tsv"
