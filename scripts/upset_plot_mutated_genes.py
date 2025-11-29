#!/usr/bin/env python3

import itertools
import matplotlib.pyplot as plt

# --- Define your gene sets per sample ---

WT_UN = ["MUC3A","MUC5AC","LILRB3","ZNF280A","CCDC187",
         "ZNF705G","NBPF26","HRNR","ANKRD36","GOLGA6L9"]

WT_APH = ["MUC3A","MUC5AC","LILRB3","ZNF280A","CCDC187",
          "ZNF705G","NBPF26","HRNR","ANKRD36","FCGBP"]

B4_UN = ["USP8","PABPC1","KMT2C","HRNR","AHNAK2","AIFM2",
         "AGAP4","ADAMTS7","ABCC8","ABCA1"]

B4_APH = ["ZFP36L1","USP8","PCDHB2","NFS1","MUC5AC",
          "MUC21","CEP170","CCT8","ABCC8","ABCA5"]

sets = {
    "WT_UN": set(WT_UN),
    "WT_APH": set(WT_APH),
    "B4_UN": set(B4_UN),
    "B4_APH": set(B4_APH),
}

labels = list(sets.keys())

# --- Compute all non-empty intersections ---

rows = []
for r in range(1, len(labels) + 1):
    for combo in itertools.combinations(labels, r):
        inter = set.intersection(*(sets[c] for c in combo))
        if len(inter) > 0:
            rows.append({"groups": combo, "size": len(inter)})

# sort by intersection size (largest first)
rows = sorted(rows, key=lambda x: x["size"], reverse=True)

sizes  = [r["size"] for r in rows]
combos = [r["groups"] for r in rows]
x      = range(len(combos))

# --- Make UpSet-style figure (bars + dot matrix) ---

fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.4)

# Top: bar chart
ax_bar = fig.add_subplot(gs[0, 0])
ax_bar.bar(x, sizes)
ax_bar.set_ylabel("Number of shared genes")
ax_bar.set_xticks([])
ax_bar.set_title("UpSet-style Intersection Plot of Most Mutated Genes")

# Bottom: dot matrix
ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_bar)

row_positions = list(range(len(labels)))

for col_idx, combo in enumerate(combos):
    present_rows = []
    for row_idx, label in enumerate(labels):
        if label in combo:
            ax_mat.scatter(col_idx, row_idx, s=60)
            present_rows.append(row_idx)
    if len(present_rows) > 1:
        ax_mat.plot([col_idx] * len(present_rows), present_rows)

ax_mat.set_yticks(row_positions)
ax_mat.set_yticklabels(labels)
ax_mat.set_xticks(list(x))
ax_mat.set_xticklabels(
    ["\n".join(c) for c in combos],
    rotation=0,
    ha="center"
)
ax_mat.set_xlabel("Sample combinations")

plt.tight_layout()
plt.subplots_adjust(bottom=0.20)
plt.savefig("upset_most_mutated_genes.png", dpi=300)
plt.savefig("upset_most_mutated_genes.pdf")  # vector for Illustrator
plt.show()


