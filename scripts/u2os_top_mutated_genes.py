#!/usr/bin/env python3
"""
Top-30 mutated genes per U2OS sample from VEP non_common TSVs.

Filtering: PICK==1 (one annotation per variant) + protein_coding BIOTYPE.
Counting: all variants per gene symbol (intronic included; large/CFS genes dominate).

Outputs (all to results/u2os/):
  top30_{sample}.csv          — ranked gene counts for each sample
  top10_combined_u2os.csv     — side-by-side top-10 table (Chapter 6 Section 6.4)
  heatmap_top_mutated_u2os.png/.pdf
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   = Path("data/u2os")
RESULT_DIR = Path("results/u2os")
FIG_DIR    = Path("figures/u2os")
RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    "WT-U": DATA_DIR / "WT-U_vep_anno_non_common.tsv.gz",
    "WT-A": DATA_DIR / "WT-A_vep_anno_non_common.tsv.gz",
    "G3-U": DATA_DIR / "G3-U_vep_anno_non_common.tsv.gz",
    "G3-A": DATA_DIR / "G3-A_vep_anno_non_common.tsv.gz",
}

TOP_N       = 30   # genes per-sample table
HEATMAP_N   = 10   # genes in heatmap / combined table


# ── per-sample counts ──────────────────────────────────────────────────────────
print("=== U2OS Top Mutated Genes ===\n")

counts = {}   # sample -> Series (SYMBOL -> n_variants), sorted descending

for sample, path in SAMPLES.items():
    print(f"Loading {sample} ...", flush=True)

    df = pd.read_csv(
        path, sep="\t",
        usecols=["PICK", "BIOTYPE", "SYMBOL", "IMPACT"],
        low_memory=False,
    )

    total_rows = len(df)
    pick1      = df[df["PICK"] == 1]
    coding     = pick1[pick1["BIOTYPE"] == "protein_coding"]

    print(f"  rows in: {total_rows:,}  →  PICK1: {len(pick1):,}  →  protein_coding: {len(coding):,}")

    per_gene = coding.groupby("SYMBOL").size().sort_values(ascending=False)
    counts[sample] = per_gene

    top30 = per_gene.head(TOP_N).reset_index()
    top30.columns = ["Gene", "n_variants"]
    top30.insert(0, "Rank", range(1, len(top30) + 1))
    out_csv = RESULT_DIR / f"top{TOP_N}_{sample}.csv"
    top30.to_csv(out_csv, index=False)

    print(f"  top-1 gene: {per_gene.index[0]} ({per_gene.iloc[0]:,} variants)")
    print(f"  saved → {out_csv}\n")


# ── combined top-10 table ──────────────────────────────────────────────────────
print("Building combined top-10 table ...")

frames = []
for sample in SAMPLES:
    top10 = counts[sample].head(HEATMAP_N).reset_index()
    top10.columns = ["Gene", f"{sample}_n"]
    top10.insert(0, "Rank", range(1, HEATMAP_N + 1))
    frames.append(top10.set_index("Rank"))

combined = pd.concat(frames, axis=1)
out_combined = RESULT_DIR / "top10_combined_u2os.csv"
combined.to_csv(out_combined)
print(f"  saved → {out_combined}\n")

# pretty print for inspection
print(combined.to_string())
print()


# ── heatmap (presence/absence across samples) ──────────────────────────────────
print("Generating heatmap ...")

sample_names = list(SAMPLES.keys())
gene_union   = []
seen         = set()
for sample in sample_names:
    for gene in counts[sample].head(HEATMAP_N).index:
        if gene not in seen:
            gene_union.append(gene)
            seen.add(gene)

n_genes   = len(gene_union)
n_samples = len(sample_names)

matrix = np.zeros((n_genes, n_samples))
for j, sample in enumerate(sample_names):
    top_set = set(counts[sample].head(HEATMAP_N).index)
    for i, gene in enumerate(gene_union):
        matrix[i, j] = 1 if gene in top_set else 0

# CFS-associated genes (known fragile-site spanning large genes)
cfs_genes = {
    "FHIT", "WWOX", "PRKN", "MACROD2", "CNTNAP2", "PTPRD",
    "LRP1B", "CSMD1", "PTPRN2", "RBFOX1", "DLGAP2", "MAGI2",
    "DPP10", "PDE4D", "CTNNA3", "ROBO2", "NRXN3",
}
is_cfs = np.array([1 if g in cfs_genes else 0 for g in gene_union]).reshape(-1, 1)

fig, (ax_main, ax_cfs) = plt.subplots(
    1, 2, figsize=(7, max(6, n_genes * 0.45)),
    gridspec_kw={"width_ratios": [4, 0.4]},
)

ax_main.imshow(matrix, aspect="auto", cmap="Greys", vmin=0, vmax=1)
ax_main.set_xticks(range(n_samples))
ax_main.set_xticklabels(sample_names, rotation=45, ha="right", fontsize=10)
ax_main.set_yticks(range(n_genes))
ax_main.set_yticklabels(gene_union, fontsize=9)
ax_main.set_title(f"Top-{HEATMAP_N} most mutated genes\n(U2OS; PICK1 protein-coding)", fontsize=11)
ax_main.set_xlabel("Sample")

for i, val in enumerate(is_cfs.flatten()):
    fc = "#333333" if val else "white"
    ax_cfs.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, edgecolor="black", facecolor=fc))
ax_cfs.set_ylim(n_genes - 0.5, -0.5)
ax_cfs.set_xlim(0, 1)
ax_cfs.axis("off")
ax_cfs.text(0.5, -0.8, "CFS", ha="center", va="top", fontsize=8)

plt.tight_layout()

for suffix in ("png", "pdf"):
    out_fig = FIG_DIR / f"heatmap_top_mutated_u2os.{suffix}"
    plt.savefig(out_fig, dpi=300, bbox_inches="tight")
    print(f"  saved → {out_fig}")

plt.close()

print("\nDone.")
