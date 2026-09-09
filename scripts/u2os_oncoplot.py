#!/usr/bin/env python3
"""
Two-panel oncoplot for U2OS coding-variant analysis.

Panel A: top-20 most-mutated coding genes (background landscape).
Panel B: genes with G3-enriched coding burden (Δcoding ≥ 3).

Outputs:
  results/u2os/top20_coding_genes.csv         (Panel A data)
  results/u2os/g3_enriched_coding_genes.csv   (Panel B data + CFS flag)
  figures/u2os/oncoplot_panel_AB_u2os.{png,pdf}
  figures/u2os/oncoplot_top_mutated_u2os.{png,pdf}  (Panel A standalone, retained)

Run from repo root:
    python3 scripts/u2os_oncoplot.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from pathlib import Path

DATA_DIR   = Path("data/u2os")
CFS_DIR    = Path("data/CFS_fixed")
RESULT_DIR = Path("results/u2os")
FIG_DIR    = Path("figures/u2os")
FIG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["WT-U", "WT-A", "G3-U", "G3-A"]
WT_SAMPLES = ["WT-U", "WT-A"]
G3_SAMPLES = ["G3-U", "G3-A"]

# ── consequence sets ───────────────────────────────────────────────────────────
CODING_CSQ = {
    "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
    "start_lost", "splice_donor_variant", "splice_acceptor_variant",
    "inframe_insertion", "inframe_deletion", "protein_altering_variant",
    "translation_start_site",
}

SEVERITY = {
    "protein_altering_variant":  1,
    "missense_variant":          2,
    "inframe_insertion":         3,
    "inframe_deletion":          3,
    "start_lost":                4,
    "stop_lost":                 4,
    "translation_start_site":    4,
    "splice_donor_variant":      5,
    "splice_acceptor_variant":   5,
    "frameshift_variant":        6,
    "stop_gained":               7,
}

CLASS_ORDER = [
    "Nonsense_Mutation", "Frameshift_Del", "Frameshift_Ins",
    "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
    "Nonstop_Mutation", "Missense_Mutation",
    "Protein_Altering", "Multi_Hit", "No_Mutation",
]

PALETTE = {
    "Missense_Mutation":  "#3F60AC",
    "Nonsense_Mutation":  "#77C576",
    "Frameshift_Del":     "#9ECAE1",
    "Frameshift_Ins":     "#F4A582",
    "Splice_Site":        "#F97F0F",
    "In_Frame_Del":       "#DBA901",
    "In_Frame_Ins":       "#B15928",
    "Nonstop_Mutation":   "#D9534F",
    "Protein_Altering":   "#984EA3",
    "Multi_Hit":          "#2DC02D",
    "No_Mutation":        "#E8E8E8",
}

SAMPLE_COLORS = ["#2166AC", "#92C5DE", "#D6604D", "#F4A582"]


def dominant_csq(csq_str, ref, alt):
    parts = set(str(csq_str).split("&"))
    best, best_sev = None, 0
    for c in parts:
        if c in SEVERITY and SEVERITY[c] > best_sev:
            best_sev, best = SEVERITY[c], c
    if best is None:
        return None
    if best == "missense_variant":          return "Missense_Mutation"
    if best == "stop_gained":               return "Nonsense_Mutation"
    if best == "frameshift_variant":
        return "Frameshift_Del" if len(str(alt)) < len(str(ref)) else "Frameshift_Ins"
    if best in ("splice_donor_variant",
                "splice_acceptor_variant"):  return "Splice_Site"
    if best == "inframe_deletion":           return "In_Frame_Del"
    if best == "inframe_insertion":          return "In_Frame_Ins"
    if best in ("stop_lost", "start_lost",
                "translation_start_site"):   return "Nonstop_Mutation"
    if best == "protein_altering_variant":   return "Protein_Altering"
    return None


def cell_class(sub):
    classes = sub["csq_class"].unique()
    if len(classes) == 0:   return "No_Mutation", 0
    if len(classes) == 1:   return classes[0], len(sub)
    return "Multi_Hit", len(sub)


# ── load flagged genes (PAR / HLA) ────────────────────────────────────────────
flagged_genes = set()
for s in SAMPLES:
    t = pd.read_csv(RESULT_DIR / f"top30_{s}.csv")
    if "flag" in t.columns:
        flagged_genes |= set(t.loc[t["flag"].notna() & (t["flag"] != ""), "Gene"])
print(f"Flagged genes excluded from selection: {sorted(flagged_genes)}\n")


# ── read TSVs ─────────────────────────────────────────────────────────────────
print("=== Reading TSVs (PICK1, protein_coding, coding CSQ, excl. chrM) ===\n")

all_variants = []

for s in SAMPLES:
    path = DATA_DIR / f"{s}_vep_anno_non_common.tsv.gz"
    print(f"Loading {s} ...", flush=True)

    reader = pd.read_csv(
        path, sep="\t",
        usecols=["CHROM", "POS", "REF", "ALT",
                 "PICK", "BIOTYPE", "SYMBOL", "Consequence"],
        low_memory=False,
        chunksize=200_000,
    )

    sample_rows = []
    for chunk in reader:
        sub = chunk[
            (chunk["CHROM"] != "chrM") &                          # exclude mitochondrial
            (chunk["PICK"] == 1) &
            (chunk["BIOTYPE"] == "protein_coding") &
            (chunk["Consequence"].apply(
                lambda c: bool(set(str(c).split("&")) & CODING_CSQ)
            ))
        ].copy()
        if len(sub):
            sample_rows.append(sub)

    df_s = pd.concat(sample_rows, ignore_index=True)
    df_s = df_s.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT", "SYMBOL"])
    df_s["csq_class"] = df_s.apply(
        lambda r: dominant_csq(r["Consequence"], r["REF"], r["ALT"]), axis=1
    )
    df_s = df_s[df_s["csq_class"].notna()].copy()
    df_s["Sample"] = s

    print(f"  coding variants (deduped, excl. chrM): {len(df_s):,}")
    print(f"  {df_s['csq_class'].value_counts().to_dict()}\n")

    all_variants.append(
        df_s[["Sample", "CHROM", "POS", "REF", "ALT", "SYMBOL", "csq_class"]]
    )

variants = pd.concat(all_variants, ignore_index=True)
print(f"Total coding variants (all samples): {len(variants):,}\n")


# ── Panel A: top-20 genes ─────────────────────────────────────────────────────
gene_totals = (
    variants[~variants["SYMBOL"].isin(flagged_genes)]
    .groupby("SYMBOL").size()
    .sort_values(ascending=False)
)
top20_genes = gene_totals.head(20).index.tolist()

print("Panel A — top-20 genes (chrM + PAR/HLA excluded):")
for i, g in enumerate(top20_genes, 1):
    print(f"  {i:2d}. {g:20s}  {gene_totals[g]}")
print()

# Build per-(gene, sample) records for Panel A
records_A = []
for gene in top20_genes:
    for sample in SAMPLES:
        sub = variants[(variants["SYMBOL"] == gene) & (variants["Sample"] == sample)]
        cls, n = cell_class(sub) if len(sub) else ("No_Mutation", 0)
        cc = sub["csq_class"].value_counts().to_dict()
        records_A.append({
            "Gene": gene, "Sample": sample, "n_variants": n,
            "Classification": cls,
            "Missense":       cc.get("Missense_Mutation", 0),
            "Nonsense":       cc.get("Nonsense_Mutation", 0),
            "Frameshift_Del": cc.get("Frameshift_Del", 0),
            "Frameshift_Ins": cc.get("Frameshift_Ins", 0),
            "Splice_Site":    cc.get("Splice_Site", 0),
            "In_Frame_Del":   cc.get("In_Frame_Del", 0),
            "In_Frame_Ins":   cc.get("In_Frame_Ins", 0),
            "Nonstop":        cc.get("Nonstop_Mutation", 0),
            "Multi_Hit":      1 if cls == "Multi_Hit" else 0,
        })

df_A = pd.DataFrame(records_A)
df_A.to_csv(RESULT_DIR / "top20_coding_genes.csv", index=False)
print(f"Saved → results/u2os/top20_coding_genes.csv\n")


# ── Panel B: G3-enriched genes (Δ ≥ 3) ───────────────────────────────────────
print("=== Panel B: computing G3 enrichment across all genes ===\n")

# Per-gene per-sample count for every non-flagged gene
counts_wide = (
    variants[~variants["SYMBOL"].isin(flagged_genes)]
    .groupby(["SYMBOL", "Sample"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=SAMPLES, fill_value=0)
)

counts_wide["mean_WT"] = counts_wide[WT_SAMPLES].mean(axis=1)
counts_wide["mean_G3"] = counts_wide[G3_SAMPLES].mean(axis=1)
counts_wide["delta"]   = counts_wide["mean_G3"] - counts_wide["mean_WT"]

panel_B_genes_df = (
    counts_wide[counts_wide["delta"] >= 3]
    .sort_values("delta", ascending=False)
    .copy()
)

# Consistent: ≥3 variants in BOTH G3 replicates
panel_B_genes_df["g3_flag"] = panel_B_genes_df.apply(
    lambda r: "consistent"
    if (r["G3-U"] >= 3 and r["G3-A"] >= 3) else "single-replicate",
    axis=1,
)

panel_B_genes = panel_B_genes_df.index.tolist()
n_B = len(panel_B_genes)
print(f"Panel B — {n_B} genes with Δ ≥ 3:")
print(panel_B_genes_df[["WT-U","WT-A","G3-U","G3-A","mean_WT","mean_G3","delta","g3_flag"]]
      .to_string())
print()

# Build per-(gene, sample) records for Panel B classification
records_B = []
for gene in panel_B_genes:
    for sample in SAMPLES:
        sub = variants[(variants["SYMBOL"] == gene) & (variants["Sample"] == sample)]
        cls, n = cell_class(sub) if len(sub) else ("No_Mutation", 0)
        cc = sub["csq_class"].value_counts().to_dict()
        records_B.append({
            "Gene": gene, "Sample": sample, "n_variants": n,
            "Classification": cls,
            "Missense":       cc.get("Missense_Mutation", 0),
            "Nonsense":       cc.get("Nonsense_Mutation", 0),
            "Frameshift_Del": cc.get("Frameshift_Del", 0),
            "Frameshift_Ins": cc.get("Frameshift_Ins", 0),
            "Splice_Site":    cc.get("Splice_Site", 0),
            "In_Frame_Del":   cc.get("In_Frame_Del", 0),
            "In_Frame_Ins":   cc.get("In_Frame_Ins", 0),
            "Nonstop":        cc.get("Nonstop_Mutation", 0),
            "Multi_Hit":      1 if cls == "Multi_Hit" else 0,
        })

df_B_cls = pd.DataFrame(records_B)


# ── CFS overlap for Panel B genes ─────────────────────────────────────────────
print("=== CFS overlap for Panel B genes ===\n")

cfs_df = pd.concat(
    [pd.read_csv(b, sep="\t", header=None,
                 names=["chr","start","end","id","score","strand"])
     for b in CFS_DIR.glob("*.bed")],
    ignore_index=True,
)

# Approximate gene coords from min/max coding-variant positions
gene_coords_B = (
    variants[variants["SYMBOL"].isin(panel_B_genes)]
    .groupby("SYMBOL")
    .agg(chrom=("CHROM", "first"), start=("POS", "min"), end=("POS", "max"))
)

def cfs_hits(gene):
    if gene not in gene_coords_B.index:
        return ""
    r   = gene_coords_B.loc[gene]
    sub = cfs_df[cfs_df["chr"] == r["chrom"]]
    ov  = sub[(r["start"] <= sub["end"]) & (r["end"] >= sub["start"])]
    return ",".join(ov["id"].tolist())

panel_B_genes_df["CFS_overlap"] = [cfs_hits(g) for g in panel_B_genes]

# Merge delta table with g3_flag and CFS into the output CSV
out_B = panel_B_genes_df[
    ["WT-U","WT-A","G3-U","G3-A","mean_WT","mean_G3","delta","g3_flag","CFS_overlap"]
].copy()
out_B.index.name = "Gene"
out_B.to_csv(RESULT_DIR / "g3_enriched_coding_genes.csv")
print(f"Saved → results/u2os/g3_enriched_coding_genes.csv\n")

print("Panel B CFS summary:")
for gene in panel_B_genes:
    row  = panel_B_genes_df.loc[gene]
    cfs  = row["CFS_overlap"] if row["CFS_overlap"] else "—"
    print(f"  {gene:20s}  Δ={row['delta']:5.1f}  {row['g3_flag']:18s}  CFS: {cfs}")
print()


# ══════════════════════════════════════════════════════════════════════════════
# Figure helpers
# ══════════════════════════════════════════════════════════════════════════════

def build_cls_matrix(df_cls, genes, samples=SAMPLES):
    mat = pd.DataFrame(index=samples, columns=genes, dtype=object)
    for _, row in df_cls.iterrows():
        mat.loc[row["Sample"], row["Gene"]] = row["Classification"]
    return mat.fillna("No_Mutation")

def draw_oncoplot(ax_bar, ax_heat, df_cls, genes, title, bar_max=None):
    """Draw bar chart (ax_bar) + heatmap (ax_heat) for a given gene list."""
    n_genes   = len(genes)
    n_samples = len(SAMPLES)

    mat_cls = build_cls_matrix(df_cls, genes)

    # counts per (gene, sample) for the bar
    cnts = (
        df_cls[df_cls["Classification"] != "No_Mutation"]
        .groupby(["Gene","Sample"])["n_variants"].sum()
        .unstack(fill_value=0)
        .reindex(index=genes, columns=SAMPLES, fill_value=0)
    )

    # ── bar ──────────────────────────────────────────────────────────────────
    bottom = np.zeros(n_genes)
    for si, (s, sc) in enumerate(zip(SAMPLES, SAMPLE_COLORS)):
        vals = cnts[s].values.astype(float)
        ax_bar.bar(range(n_genes), vals, bottom=bottom,
                   color=sc, alpha=0.85, width=0.7, label=s)
        bottom += vals

    ax_bar.set_xlim(-0.5, n_genes - 0.5)
    ax_bar.set_xticks([])
    ax_bar.set_ylabel("Coding\nvariants", fontsize=7.5)
    if bar_max:
        ax_bar.set_ylim(0, bar_max)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.tick_params(labelsize=7)
    ax_bar.set_title(title, fontsize=9.5, loc="left", pad=4, fontweight="bold")
    ax_bar.legend(title="Sample", fontsize=6, title_fontsize=6.5,
                  loc="upper right", frameon=False, ncol=2)

    # ── heatmap ───────────────────────────────────────────────────────────────
    for gi, gene in enumerate(genes):
        for si, sample in enumerate(SAMPLES):
            cls   = mat_cls.loc[sample, gene]
            color = PALETTE.get(cls, PALETTE["No_Mutation"])
            ax_heat.add_patch(
                plt.Rectangle((gi - 0.45, si - 0.45), 0.9, 0.9,
                               facecolor=color, edgecolor="white", linewidth=0.5)
            )
            if cls == "Multi_Hit":
                ax_heat.plot(gi, si, "w.", markersize=5, zorder=10)

    ax_heat.set_xlim(-0.5, n_genes - 0.5)
    ax_heat.set_ylim(-0.5, n_samples - 0.5)
    ax_heat.set_xticks(range(n_genes))
    ax_heat.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticks(range(n_samples))
    ax_heat.set_yticklabels(SAMPLES, fontsize=8.5)
    ax_heat.axhline(1.5, color="#999", lw=0.7, ls="--")
    ax_heat.text(-0.6, 0.5, "WT", ha="right", va="center",
                 fontsize=7.5, color="#444", style="italic")
    ax_heat.text(-0.6, 2.5, "G3", ha="right", va="center",
                 fontsize=7.5, color="#c0392b", style="italic", fontweight="bold")
    ax_heat.set_frame_on(False)
    ax_heat.tick_params(length=0)

    # Δ annotations above each Panel B gene column (if delta info available)
    return mat_cls


def draw_legend(ax, all_cls_values):
    ax.axis("off")
    present = [c for c in CLASS_ORDER if c in all_cls_values]
    patches = [
        mpatches.Patch(facecolor=PALETTE[c], edgecolor="#999",
                       linewidth=0.4, label=c.replace("_", " "))
        for c in present
    ]
    ax.legend(handles=patches, title="Classification",
              title_fontsize=8, fontsize=7.5,
              loc="upper left", frameon=True, framealpha=0.9,
              edgecolor="#ccc", borderpad=0.9)


# ══════════════════════════════════════════════════════════════════════════════
# Two-panel figure  (stacked: A above B, shared legend column)
# ══════════════════════════════════════════════════════════════════════════════
print("=== Generating two-panel figure ===\n")

# Proportional widths: Panel A always 20 genes wide; Panel B has n_B genes.
# We normalise so Panel A occupies the full column width.
# For layout, we fix gene column widths at 0.55 inches each.
GENE_W    = 0.52    # inches per gene
BAR_H     = 0.9     # inches for bar panel
HEAT_H    = 1.5     # inches for heatmap (4 samples)
GAP_H     = 0.6     # gap between panels
LEG_W     = 2.0     # legend column width
PAD_L     = 1.1     # left margin (for sample labels)
PAD_R     = 0.2

fig_w = PAD_L + max(20, n_B) * GENE_W + LEG_W + PAD_R
fig_h = BAR_H + HEAT_H + GAP_H + BAR_H + HEAT_H + 0.8   # 0.8 top/bottom margin

fig = plt.figure(figsize=(fig_w, fig_h))

content_w = max(20, n_B) * GENE_W
leg_frac  = LEG_W / fig_w
con_frac  = 1.0 - leg_frac

# GridSpec: 5 rows (barA, heatA, gap, barB, heatB) × 2 cols (content, legend)
gs = gridspec.GridSpec(
    5, 2,
    figure=fig,
    height_ratios=[BAR_H, HEAT_H, GAP_H, BAR_H, HEAT_H],
    width_ratios=[con_frac, leg_frac],
    hspace=0.0,
    wspace=0.03,
    left=PAD_L / fig_w,
    right=1.0 - PAD_R / fig_w,
    top=0.97,
    bottom=0.08,
)

ax_A_bar  = fig.add_subplot(gs[0, 0])
ax_A_heat = fig.add_subplot(gs[1, 0])
ax_B_bar  = fig.add_subplot(gs[3, 0])
ax_B_heat = fig.add_subplot(gs[4, 0])
ax_leg    = fig.add_subplot(gs[:, 1])

# ── draw Panel A ──────────────────────────────────────────────────────────────
mat_A = draw_oncoplot(
    ax_A_bar, ax_A_heat, df_A, top20_genes,
    title="A.  Most-mutated coding genes (background landscape)",
)

# ── draw Panel B ──────────────────────────────────────────────────────────────
# Add Δ value to gene labels
delta_labels = [
    f"{g}\n(Δ={panel_B_genes_df.loc[g,'delta']:.0f})" for g in panel_B_genes
]
# Substitute in ax after draw so we can modify xticklabels
mat_B = draw_oncoplot(
    ax_B_bar, ax_B_heat, df_B_cls, panel_B_genes,
    title=f"B.  Coding genes enriched in ZFP36L1−/−  (Δcoding ≥ 3,  n = {n_B})",
)

# Override x-tick labels to show Δ and consistent flag
xtick_labels_B = []
for g in panel_B_genes:
    row   = panel_B_genes_df.loc[g]
    flag  = "●" if row["g3_flag"] == "consistent" else "○"   # ● = consistent
    label = f"{g}\nΔ={row['delta']:.0f} {flag}"
    xtick_labels_B.append(label)

ax_B_heat.set_xticklabels(xtick_labels_B, rotation=45, ha="right", fontsize=7.5)

# If Panel B has fewer genes than 20, pad the x-axis so cells are square-ish
# (already handled by the fixed GENE_W in figure sizing)
# Shrink Panel B x-axis to match actual gene count (leave rest blank)
ax_B_bar.set_xlim(-0.5, n_B - 0.5)
ax_B_heat.set_xlim(-0.5, n_B - 0.5)
ax_A_bar.set_xlim(-0.5, 19.5)
ax_A_heat.set_xlim(-0.5, 19.5)

# ── legend (shared) ───────────────────────────────────────────────────────────
all_cls = set(mat_A.values.flatten()) | set(mat_B.values.flatten())
draw_legend(ax_leg, all_cls)

# Consistent/single-replicate key below the legend
ax_leg.text(0.05, 0.35,
            "● consistent (≥3 in both G3)\n○ single-replicate",
            transform=ax_leg.transAxes,
            fontsize=7, va="top", color="#333")

for suffix in ("png", "pdf"):
    out = FIG_DIR / f"oncoplot_panel_AB_u2os.{suffix}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved → {out}")

plt.close()

# ── also save standalone Panel A (replaces old file) ─────────────────────────
fig2, axes2 = plt.subplots(2, 2,
    figsize=(14, 5.5),
    gridspec_kw={"height_ratios":[1,3], "width_ratios":[20,2], "hspace":0.06, "wspace":0.03}
)
draw_oncoplot(axes2[0,0], axes2[1,0], df_A, top20_genes,
              title="Top-20 mutated coding genes — U2OS")
draw_legend(axes2[0,1], set(build_cls_matrix(df_A, top20_genes).values.flatten()))
axes2[1,1].axis("off")
for suffix in ("png","pdf"):
    plt.savefig(FIG_DIR / f"oncoplot_top_mutated_u2os.{suffix}", dpi=300, bbox_inches="tight")
plt.close()
print(f"Updated standalone oncoplot_top_mutated_u2os.{{png,pdf}}\n")

print("Done.")
