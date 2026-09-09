#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

# -------------------------
# Data: gene lists per sample
# -------------------------
WT_UN = ["MUC3A","MUC5AC","LILRB3","ZNF280A","CCDC187",
         "ZNF705G","NBPF26","HRNR","ANKRD36","GOLGA6L9"]

WT_APH = ["MUC3A","MUC5AC","LILRB3","ZNF280A","CCDC187",
          "ZNF705G","NBPF26","HRNR","ANKRD36","FCGBP"]

B4_UN  = ["USP8","PABPC1","KMT2C","HRNR","AHNAK2","AIFM2",
          "AGAP4","ADAMTS7","ABCC8","ABCA1"]

B4_APH = ["ZFP36L1","USP8","PCDHB2","NFS1","MUC5AC",
          "MUC21","CEP170","CCT8","ABCC8","ABCA5"]

samples = ["WT_UN", "WT_APH", "B4_UN", "B4_APH"]
gene_sets = {"WT_UN": WT_UN, "WT_APH": WT_APH, "B4_UN": B4_UN, "B4_APH": B4_APH}

# -------------------------
# Build gene list & matrix
# -------------------------
genes = sorted(list(set(WT_UN + WT_APH + B4_UN + B4_APH)))

matrix = np.zeros((len(genes), len(samples)))
for i, gene in enumerate(genes):
    for j, sample in enumerate(samples):
        matrix[i, j] = 1 if gene in gene_sets[sample] else 0

# -------------------------
# Gene functional categories
# -------------------------
category_map = {
    "MUC3A": "Mucin",
    "MUC5AC": "Mucin",
    "MUC21": "Mucin",
    "ZNF280A": "ZNF",
    "ZNF705G": "ZNF",
    "ABCC8": "Transporter",
    "ABCA1": "Transporter",
    "ABCA5": "Transporter",
    "USP8": "Other",
    "PABPC1": "Other",
    "KMT2C": "Chromatin",
    "HRNR": "Scaffold",
    "AIFM2": "Apoptosis",
    "AGAP4": "Trafficking",
    "ADAMTS7": "Protease",
    "GOLGA6L9": "Other",
    "NBPF26": "Other",
    "CCDC187": "Other",
    "FCGBP": "Other",
    "NFS1": "Other",
    "PCDHB2": "Adhesion",
    "CEP170": "Centrosome",
    "CCT8": "Chaperone",
    "ZFP36L1": "RNA_binding",
}

categories = [category_map.get(g, "Other") for g in genes]
cat_types = sorted(set(categories))
cat_colors = {
    "Mucin": "red",
    "ZNF": "blue",
    "Transporter": "green",
    "Chromatin": "purple",
    "Scaffold": "orange",
    "Apoptosis": "brown",
    "Trafficking": "gold",
    "Protease": "pink",
    "Adhesion": "cyan",
    "Centrosome": "olive",
    "Chaperone": "magenta",
    "RNA_binding": "teal",
    "Other": "lightgrey",
}
cat_color_list = [cat_colors[c] for c in categories]

# -------------------------
# Fragile-site membership
# -------------------------
cfs_genes = {
    "ABCC8","AGAP4","AIFM2","ANKRD36","FCGBP",
    "HRNR","KMT2C","LILRB3","MUC21","MUC3A","ZFP36L1"
}
is_cfs = np.array([1 if g in cfs_genes else 0 for g in genes]).reshape(-1,1)

# -------------------------
# Plot: main heatmap + category bar + CFS bar
# -------------------------
fig, (ax_main, ax_cat, ax_cfs) = plt.subplots(
    1, 3, figsize=(10, 14),
    gridspec_kw={"width_ratios": [4, 0.4, 0.4]}
)

# Main presence/absence heatmap
ax_main.imshow(matrix, aspect="auto", cmap="Greys")
ax_main.set_xticks(range(len(samples)))
ax_main.set_xticklabels(samples, rotation=45, ha="right")
ax_main.set_yticks(range(len(genes)))
ax_main.set_yticklabels(genes)
ax_main.set_title("Most Mutated Genes Across Samples")
ax_main.set_xlabel("Samples")
ax_main.set_ylabel("Genes")

# Category annotation bar
for i, color in enumerate(cat_color_list):
    ax_cat.add_patch(plt.Rectangle((0, i-0.5), 1, 1, color=color))
ax_cat.set_ylim(len(genes)-0.5, -0.5)
ax_cat.set_xlim(0, 1)
ax_cat.axis("off")

used_cats = sorted(set(categories))
legend_handles = [
    plt.Line2D([0],[0], marker='s', linestyle='None',
               markersize=8, markerfacecolor=cat_colors[c],
               label=c)
    for c in used_cats
]
ax_cat.legend(handles=legend_handles, title="Category",
              loc="upper left", bbox_to_anchor=(1.05, 1.0))

# CFS annotation bar
for i, val in enumerate(is_cfs.flatten()):
    color = "black" if val == 1 else "white"
    ax_cfs.add_patch(
        plt.Rectangle((0, i-0.5), 1, 1,
                      edgecolor="black", facecolor=color)
    )
ax_cfs.set_ylim(len(genes)-0.5, -0.5)
ax_cfs.set_xlim(0, 1)
ax_cfs.axis("off")
ax_cfs.text(0.5, -1, "CFS", ha="center", va="top")

plt.tight_layout()
plt.savefig("heatmap_mutated_genes_annotated.png", dpi=300)
plt.savefig("heatmap_mutated_genes_annotated.pdf")
plt.show()


