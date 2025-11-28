# Replication Stress WGS Analysis

This repository contains the computational workflow, Python scripts, 
and reproducible analysis pipeline used to generate the figures and 
tables for the thesis chapters on replication stress, CRISPR-mediated 
ZFP36L1 knockout, fragile-site mapping, and whole-genome sequencing.

## Project structure

replication-stress-wgs-analysis/
│
├── scripts/ # Python analysis scripts (UpSet plot, CFS mapping, heatmaps)
├── data/ # Input datasets (e.g. gene lists, BED files)
├── results/ # Output tables and processed data
├── figures/ # Generated figures (PNG, PDF, SVG)
└── notebooks/ # Optional Jupyter notebooks


## Objectives

- Identify and compare mutation landscapes across cell lines and conditions
- Map mutated genes to common fragile sites (CFS)
- Visualise overlapping mutation patterns using UpSet plots
- Reproduce the analyses used in the PhD thesis

## Requirements

- Python 3.10+
- matplotlib
- pandas

## Usage

Each script in `scripts/` can be run independently.  
Example:

```bash
python scripts/upset_plot_mutated_genes.py

