#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import zipfile
import glob
import os

# ---------------------------------------------------------
# 1. Load mutated genes table
# ---------------------------------------------------------

genes = pd.read_excel("data/MCF7MutatedGenesLocations.xlsx", header=1)

# remove stray index column if present
if "Unnamed: 0" in genes.columns:
    genes = genes.drop(columns=["Unnamed: 0"])

genes["Gene_Start"] = pd.to_numeric(genes["Gene_Start"])
genes["Gene_End"]   = pd.to_numeric(genes["Gene_End"])
genes["CHROMOSOME"] = genes["CHROMOSOME"].astype(str)


# ---------------------------------------------------------
# 2. Load CFS sites from cfsfixed.zip
# ---------------------------------------------------------

extract_folder = "data/cfs_extracted"

if not os.path.exists(extract_folder):
    os.makedirs(extract_folder)

with zipfile.ZipFile("data/cfsfixed.zip") as z:
    z.extractall(extract_folder)

cfs_list = []
for bf in glob.glob(os.path.join(extract_folder, "*.bed")):
    df = pd.read_csv(
        bf, sep="\t", header=None,
        names=["chr", "start", "end", "id", "score", "strand"]
    )
    cfs_list.append(df)

cfs = pd.concat(cfs_list, ignore_index=True)


# ---------------------------------------------------------
# 3. Chromosome lengths for hg38
# ---------------------------------------------------------

chrom_lengths = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345,  "chr17": 83257441,  "chr18": 80373285,
    "chr19": 58617616,  "chr20": 64444167,  "chr21": 46709983,
    "chr22": 50818468,  "chrX": 156040895
}

# X chromosome placed at far right
chrom_order = list(chrom_lengths.keys())


# ---------------------------------------------------------
# 4. Simple vertical label de-overlap
# ---------------------------------------------------------

def resolve_overlaps(yvals, min_sep=3e6):
    """
    Given a list of y positions (midpoints), enforce a minimum
    vertical separation. Used per chromosome.
    """
    adjusted = yvals.copy()
    for i in range(1, len(adjusted)):
        if adjusted[i] - adjusted[i-1] < min_sep:
            adjusted[i] = adjusted[i-1] + min_sep
    return adjusted


# ---------------------------------------------------------
# 5. Plotting
# ---------------------------------------------------------

fig, ax = plt.subplots(figsize=(22, 12))

xpos = {chrom: i for i, chrom in enumerate(chrom_order)}

label_shift = 0.20      # horizontal shift of labels
font_size   = 10        # label font size
min_sep     = 3e6       # vertical spacing for de-overlap


for chrom in chrom_order:

    x = xpos[chrom]
    chr_len = chrom_lengths[chrom]

    # Draw chromosome backbone
    ax.plot([x, x], [0, chr_len], color="lightgrey", linewidth=8, zorder=1)

    # Draw CFS blocks
    csub = cfs[cfs["chr"] == chrom]
    for _, r in csub.iterrows():
        ax.add_patch(
            plt.Rectangle(
                (x - 0.05, r["start"]),
                0.10,
                r["end"] - r["start"],
                facecolor="black",
                edgecolor="black",
                zorder=3
            )
        )

    # Select mutated genes on this chromosome
    gsub = genes[genes["CHROMOSOME"] == chrom].copy()

    # One label per gene name
    gsub = gsub.drop_duplicates(subset=["GENE_NAME"], keep="first")

    # Midpoints
    gsub["mid"] = (gsub["Gene_Start"] + gsub["Gene_End"]) / 2
    gsub = gsub.sort_values("mid")

    orig_positions = list(gsub["mid"])
    adj_positions  = resolve_overlaps(orig_positions, min_sep=min_sep)

    # Draw labels + leader lines
    for (_, row), new_y in zip(gsub.iterrows(), adj_positions):

        # Leader line (midpoint → adjusted label position)
        ax.annotate(
            "",
            xy=(x, row["mid"]),
            xytext=(x + label_shift, new_y),
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            zorder=2
        )

        # Gene label
        ax.text(
            x + label_shift,
            new_y,
            row["GENE_NAME"],
            fontsize=font_size,
            va="center",
            zorder=10
        )


# ---------------------------------------------------------
# 6. Final formatting
# ---------------------------------------------------------

ax.set_xlim(-1, len(chrom_order) + 1)
ax.set_xticks(range(len(chrom_order)))
ax.set_xticklabels([c.replace("chr", "") for c in chrom_order])
ax.set_ylabel("Genomic position (bp)")
ax.set_title("Genome-wide Ideogram — Unique Gene Labels with Leader Lines", fontsize=16)

plt.tight_layout()

# ---------------------------------------------------------
# 6. Save figure to disk
# ---------------------------------------------------------

# Create output directory if needed
output_dir = "figures"
os.makedirs(output_dir, exist_ok=True)

output_path = os.path.join(output_dir, "genomewide_ideogram_mutated_genes.png")

plt.tight_layout()
plt.savefig(output_path, dpi=300, bbox_inches="tight")

print(f"Saved figure to: {output_path}")

plt.show()


