#!/usr/bin/env python3
"""
Plot a genome-wide ideogram of common fragile sites (CFSs) on true hg38
chromosome lengths, using per-chromosome BED files.

Colours:
- red / green  : adjacent, non-overlapping CFS blocks (alternate)
- black        : CFS blocks that overlap the previous block on that chromosome

Input:
    data/CFS_fixed/chr*_fragile_site.bed

Each .bed file is expected to have columns:
    chrom  start  end  name  score  strand

Outputs:
    figures/CFS_ideogram_hg38_scaled.png
    figures/CFS_ideogram_hg38_scaled.pdf
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

# Folder containing per-chromosome BED files (adjust if needed)
BED_FOLDER = "data/CFS_fixed"

# Output folder for figures
FIG_DIR = "figures"

# Chromosomes to plot (autosomes + X)
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

# hg38 chromosome lengths in base pairs
HG38_CHROM_LENGTHS = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415,
}


# ---------------------------------------------------------------------
# LOAD CFS DATA
# ---------------------------------------------------------------------

all_rows = []
for chrom in CHROM_ORDER:
    bed_path = os.path.join(BED_FOLDER, f"{chrom}_fragile_site.bed")
    if not os.path.exists(bed_path):
        continue
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    all_rows.append(df)

if not all_rows:
    raise RuntimeError(f"No .bed files found in {BED_FOLDER}")

cfs = pd.concat(all_rows, ignore_index=True)

# Convert to megabases for plotting
cfs["start_mb"] = cfs["start"] / 1e6
cfs["end_mb"] = cfs["end"] / 1e6


# ---------------------------------------------------------------------
# ASSIGN COLOURS PER CHROMOSOME (overlap logic)
# ---------------------------------------------------------------------

def assign_colors(df_chr: pd.DataFrame) -> pd.DataFrame:
    """Assign red/green/black colours for one chromosome.

    - First block: red
    - If current.start <= previous.end  -> black (overlap)
    - Else alternate red/green for adjacent non-overlapping blocks
    """
    df_chr = df_chr.sort_values("start").reset_index(drop=True)
    colors = []
    for i in range(len(df_chr)):
        if i == 0:
            colors.append("black")
        else:
            prev_end = df_chr.loc[i - 1, "end_mb"]
            curr_start = df_chr.loc[i, "start_mb"]
            if curr_start <= prev_end:
                colors.append("black")  # overlapping with previous
            else:
                colors.append("black" if colors[-1] == "black" else "black")
    df_chr["color"] = colors
    return df_chr


colored_rows = []
for chrom in CHROM_ORDER:
    sub = cfs[cfs["chrom"] == chrom]
    if len(sub) > 0:
        colored_rows.append(assign_colors(sub))

cfs_colored = pd.concat(colored_rows, ignore_index=True)


# ---------------------------------------------------------------------
# PLOT HG38-SCALED IDEOGRAM
# ---------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(12, 6))
x_positions = range(len(CHROM_ORDER))

for x, chrom in zip(x_positions, CHROM_ORDER):
    chrom_len_bp = HG38_CHROM_LENGTHS[chrom]
    chrom_len_mb = chrom_len_bp / 1e6

    # Draw chromosome backbone
    ax.plot([x, x], [0, chrom_len_mb], color="lightgrey", linewidth=4, zorder=1)

    # Draw CFS blocks on that chromosome
    sub = cfs_colored[cfs_colored["chrom"] == chrom]
    for _, row in sub.iterrows():
        ax.add_patch(
            plt.Rectangle(
                (x - 0.25, row["start_mb"]),             # x, y
                0.5,                                     # width
                row["end_mb"] - row["start_mb"],         # height
                facecolor=row["color"],
                edgecolor="black",
                linewidth=0.8,
                alpha=0.9,
                zorder=2,
            )
        )

ax.set_xticks(list(x_positions))
ax.set_xticklabels([c.replace("chr", "") for c in CHROM_ORDER])
ax.set_xlim(-0.5, len(CHROM_ORDER) - 0.5)
ax.set_ylabel("Genomic position (Mb)")
ax.set_xlabel("Chromosome")
ax.set_title("Genome-wide Common Fragile Sites on hg38 Chromosome-Length Ideogram")

plt.tight_layout()

os.makedirs(FIG_DIR, exist_ok=True)
out_png = os.path.join(FIG_DIR, "CFS_ideogram_hg38_scaled.png")
out_pdf = os.path.join(FIG_DIR, "CFS_ideogram_hg38_scaled.pdf")

plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf, bbox_inches="tight")
plt.show()

print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")

