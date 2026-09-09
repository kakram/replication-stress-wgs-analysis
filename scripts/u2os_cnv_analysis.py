#!/usr/bin/env python3
"""
CNV analysis — U2OS WGS (Section 6.5).

Inputs:  data/u2os/{WT-U,WT-A,G3-U,G3-A}_cnv.cns   (CNVkit segment files)
Outputs:
  results/u2os/cnv_per_chromosome_summary.csv
  results/u2os/cnv_targeted_loci.csv
  results/u2os/cnv_genome_burden.csv
  figures/u2os/cnv_genomewide_u2os.{png,pdf}

Run from repo root:
    python3 scripts/u2os_cnv_analysis.py

Threshold definitions (Methods):
  log2 < -1.0               → deletion
  -1.0 ≤ log2 < -0.25      → loss
  -0.25 ≤ log2 ≤  0.25     → neutral
   0.25 < log2 ≤  1.0      → gain
   log2 >  1.0              → amplification
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

DATA_DIR   = Path("data/u2os")
RESULT_DIR = Path("results/u2os")
FIG_DIR    = Path("figures/u2os")
RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["WT-U", "WT-A", "G3-U", "G3-A"]

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# ── CNV call thresholds ────────────────────────────────────────────────────────
def call_segment(log2):
    if   log2 <  -1.00: return "deletion"
    elif log2 <  -0.25: return "loss"
    elif log2 <=  0.25: return "neutral"
    elif log2 <=  1.00: return "gain"
    else:               return "amplification"

CALL_COLORS = {
    "deletion":      "#08306B",   # dark navy
    "loss":          "#4292C6",   # medium blue
    "neutral":       "#CCCCCC",   # grey
    "gain":          "#FC8D59",   # orange-red
    "amplification": "#D73027",   # dark red
}

# For the genome-wide plot: collapse to 3 colours (user spec)
PLOT_COLOR = {
    "deletion":      "#2166AC",   # blue family
    "loss":          "#6BAED6",   # blue family
    "neutral":       "#CCCCCC",
    "gain":          "#EF3B2C",   # red family
    "amplification": "#99000D",   # red family
}


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load & concatenate CNS files
# ══════════════════════════════════════════════════════════════════════════════
print("=== Loading CNS files ===\n")

frames = []
for s in SAMPLES:
    path = DATA_DIR / f"{s}_cnv.cns"
    df   = pd.read_csv(path, sep="\t")
    df["sample"] = s
    df["length"] = df["end"] - df["start"]
    df["call"]   = df["log2"].apply(call_segment)
    frames.append(df)

segs = pd.concat(frames, ignore_index=True)

# Keep only canonical chromosomes
segs = segs[segs["chromosome"].isin(CHROM_ORDER)].copy()

print(f"Total segments: {len(segs):,}  across {segs['sample'].nunique()} samples")
for s in SAMPLES:
    sub = segs[segs["sample"] == s]
    print(f"  {s}: {len(sub)} segments  "
          f"log2 [{sub['log2'].min():.2f}, {sub['log2'].max():.2f}]")
print()


# ══════════════════════════════════════════════════════════════════════════════
# 2. Per-chromosome summary
# ══════════════════════════════════════════════════════════════════════════════
print("=== Per-chromosome summary ===\n")

chr_rows = []
for chrom in CHROM_ORDER:
    sub_c = segs[segs["chromosome"] == chrom]
    if sub_c.empty:
        continue
    for s in SAMPLES:
        sub = sub_c[sub_c["sample"] == s]
        if sub.empty:
            continue
        total_len = sub["length"].sum()
        wmean_log2 = (sub["log2"] * sub["length"]).sum() / total_len

        for call in ["deletion","loss","neutral","gain","amplification"]:
            bp = sub.loc[sub["call"] == call, "length"].sum()
            chr_rows.append({
                "chromosome": chrom,
                "sample":     s,
                "call":       call,
                "Mb":         round(bp / 1e6, 3),
                "pct":        round(100 * bp / total_len, 2),
                "wmean_log2": round(wmean_log2, 4),
            })

chr_summary = pd.DataFrame(chr_rows)
chr_summary.to_csv(RESULT_DIR / "cnv_per_chromosome_summary.csv", index=False)
print(f"Saved → results/u2os/cnv_per_chromosome_summary.csv\n")


# ══════════════════════════════════════════════════════════════════════════════
# 3. Genome-wide CNV burden per sample
# ══════════════════════════════════════════════════════════════════════════════
print("=== Genome-wide CNV burden ===\n")

burden_rows = []
for s in SAMPLES:
    sub = segs[segs["sample"] == s]
    total_mb = sub["length"].sum() / 1e6
    row = {"sample": s, "total_Mb": round(total_mb, 1)}
    for call in ["deletion","loss","neutral","gain","amplification"]:
        mb = sub.loc[sub["call"] == call, "length"].sum() / 1e6
        row[f"{call}_Mb"]  = round(mb, 1)
        row[f"{call}_pct"] = round(100 * mb / total_mb, 1)
    burden_rows.append(row)

burden_df = pd.DataFrame(burden_rows)
burden_df.to_csv(RESULT_DIR / "cnv_genome_burden.csv", index=False)

print(burden_df[["sample","deletion_Mb","loss_Mb","neutral_Mb",
                  "gain_Mb","amplification_Mb"]].to_string(index=False))
print()


# ══════════════════════════════════════════════════════════════════════════════
# 4. Genome-wide CNV figure
# ══════════════════════════════════════════════════════════════════════════════
print("=== Generating genome-wide figure ===\n")

# Chromosome sizes and offsets (from data)
chr_sizes = (
    segs.groupby("chromosome")["end"].max()
    .reindex(CHROM_ORDER).dropna().astype(int)
)
chr_offsets = {}
offset = 0
for chrom in chr_sizes.index:
    chr_offsets[chrom] = offset
    offset += int(chr_sizes[chrom])
genome_size = offset

PLOT_Y_MIN, PLOT_Y_MAX = -2.0, 2.0

fig, axes = plt.subplots(
    len(SAMPLES), 1,
    figsize=(16, 7),
    sharex=True,
    gridspec_kw={"hspace": 0.12}
)

for ax, s in zip(axes, SAMPLES):
    sub = segs[segs["sample"] == s].copy()
    sub = sub[sub["chromosome"].isin(chr_offsets)].copy()
    sub["x_start"] = sub["chromosome"].map(chr_offsets) + sub["start"]
    sub["x_end"]   = sub["chromosome"].map(chr_offsets) + sub["end"]
    sub["log2_clipped"] = sub["log2"].clip(PLOT_Y_MIN, PLOT_Y_MAX)

    for _, row in sub.iterrows():
        color = PLOT_COLOR[row["call"]]
        # Draw filled rectangle from 0 to log2 (gains up, losses down)
        y0 = min(0.0, row["log2_clipped"])
        h  = abs(row["log2_clipped"])
        ax.add_patch(plt.Rectangle(
            (row["x_start"], y0), row["x_end"] - row["x_start"], h,
            facecolor=color, edgecolor="none", alpha=0.85
        ))

    # Baseline
    ax.axhline(0, color="#666", lw=0.5, zorder=5)
    # Threshold lines
    for thresh, ls in [(-0.25, ":"), (0.25, ":"), (-1.0, "--"), (1.0, "--")]:
        ax.axhline(thresh, color="#999", lw=0.4, ls=ls, zorder=4)

    # Chromosome separators and labels
    for chrom in chr_sizes.index:
        x = chr_offsets[chrom]
        ax.axvline(x, color="#BBBBBB", lw=0.4, zorder=3)
        # Label at midpoint of chromosome
        mid = x + chr_sizes[chrom] / 2
        label = chrom.replace("chr", "")
        ax.text(mid, PLOT_Y_MAX * 0.88, label,
                ha="center", va="top", fontsize=5.5, color="#555")

    ax.set_xlim(0, genome_size)
    ax.set_ylim(PLOT_Y_MIN, PLOT_Y_MAX)
    ax.set_ylabel(s, fontsize=8.5, rotation=0, labelpad=35, va="center")
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(["-2", "-1", "0", "+1", "+2"], fontsize=6.5)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

axes[-1].set_xlabel("Genomic position (chr1 → chrY)", fontsize=8)
axes[-1].tick_params(axis="x", bottom=True, labelbottom=False)

# Legend
legend_patches = [
    mpatches.Patch(color=PLOT_COLOR["amplification"], label="Amplification (log2 > 1)"),
    mpatches.Patch(color=PLOT_COLOR["gain"],          label="Gain (0.25 < log2 ≤ 1)"),
    mpatches.Patch(color=PLOT_COLOR["neutral"],       label="Neutral"),
    mpatches.Patch(color=PLOT_COLOR["loss"],          label="Loss (−1 ≤ log2 < −0.25)"),
    mpatches.Patch(color=PLOT_COLOR["deletion"],      label="Deletion (log2 < −1)"),
]
axes[0].legend(handles=legend_patches, fontsize=6.5, loc="upper right",
               frameon=True, framealpha=0.9, edgecolor="#ccc",
               title="CNV call", title_fontsize=7)

fig.suptitle("Genome-wide CNV — U2OS (CNVkit segments, GRCh38)",
             fontsize=10, y=1.01)

for suffix in ("png", "pdf"):
    out = FIG_DIR / f"cnv_genomewide_u2os.{suffix}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved → {out}")
plt.close()


# ══════════════════════════════════════════════════════════════════════════════
# 5. Targeted locus verification
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Targeted locus verification ===\n")

LOCI = {
    "FCGBP":   ("chr19",  39_872_500,  39_966_000),
    "FAM8A1":  ("chr6",   17_602_000,  17_673_000),
    "ZNF208":  ("chr19",  21_994_000,  22_073_000),
    "FHIT":    ("chr3",   59_747_000,  61_237_000),
    "ROBO2":   ("chr3",   75_907_852,  77_651_517),
    "ZFP36L1": ("chr14",  68_792_000,  68_797_000),
}

locus_rows = []
for gene, (chrom, lstart, lend) in LOCI.items():
    for s in SAMPLES:
        sub = segs[
            (segs["sample"] == s) &
            (segs["chromosome"] == chrom) &
            (segs["start"] <= lend) &
            (segs["end"]   >= lstart)
        ]
        if sub.empty:
            locus_rows.append({
                "gene": gene, "sample": s, "locus": f"{chrom}:{lstart:,}-{lend:,}",
                "n_segments": 0,
                "log2_values": "—", "calls": "—",
                "dominant_call": "no_data", "wmean_log2": np.nan,
            })
        else:
            # Weighted mean log2 by overlap length
            overlap_len = (
                sub.apply(lambda r: min(r["end"], lend) - max(r["start"], lstart), axis=1)
            )
            wmean = (sub["log2"] * overlap_len).sum() / overlap_len.sum()
            locus_rows.append({
                "gene":          gene,
                "sample":        s,
                "locus":         f"{chrom}:{lstart:,}-{lend:,}",
                "n_segments":    len(sub),
                "log2_values":   ", ".join(f"{v:.3f}" for v in sub["log2"]),
                "calls":         ", ".join(sub["call"]),
                "dominant_call": call_segment(wmean),
                "wmean_log2":    round(wmean, 4),
            })

locus_df = pd.DataFrame(locus_rows)
locus_df.to_csv(RESULT_DIR / "cnv_targeted_loci.csv", index=False)
print(f"Saved → results/u2os/cnv_targeted_loci.csv\n")

# Pretty print
for gene in LOCI:
    chrom, lstart, lend = LOCI[gene]
    print(f"  {gene}  ({chrom}:{lstart:,}–{lend:,})")
    sub = locus_df[locus_df["gene"] == gene]
    for _, r in sub.iterrows():
        print(f"    {r['sample']:6s}  log2={r['wmean_log2']:+.3f}  call={r['dominant_call']:14s}  "
              f"({r['n_segments']} segment(s): {r['calls']})")
    print()


# ══════════════════════════════════════════════════════════════════════════════
# 6. Chromosomal redistribution: chr3, chr15, chr18
# ══════════════════════════════════════════════════════════════════════════════
print("=== Chr3 / Chr15 / Chr18 — redistribution vs CNV ===\n")

focus_chroms = ["chr3", "chr15", "chr18"]

# Weighted mean log2 per chromosome per sample
chrom_wmean = {}
for chrom in focus_chroms:
    chrom_wmean[chrom] = {}
    for s in SAMPLES:
        sub = segs[(segs["chromosome"] == chrom) & (segs["sample"] == s)]
        if sub.empty:
            chrom_wmean[chrom][s] = np.nan
            continue
        total = sub["length"].sum()
        chrom_wmean[chrom][s] = (sub["log2"] * sub["length"]).sum() / total

print(f"{'Chromosome':10s}  {'WT-U':>7s}  {'WT-A':>7s}  {'WT mean':>8s}  "
      f"{'G3-U':>7s}  {'G3-A':>7s}  {'G3 mean':>8s}  {'Δ(G3-WT)':>9s}  Verdict")
print("-" * 95)

redistrib_rows = []
for chrom in focus_chroms:
    wt_u = chrom_wmean[chrom]["WT-U"]
    wt_a = chrom_wmean[chrom]["WT-A"]
    g3_u = chrom_wmean[chrom]["G3-U"]
    g3_a = chrom_wmean[chrom]["G3-A"]
    wt_mean = np.mean([wt_u, wt_a])
    g3_mean = np.mean([g3_u, g3_a])
    delta   = g3_mean - wt_mean

    if   abs(delta) < 0.10: verdict = "no clear CNV difference"
    elif delta >  0.25:     verdict = "CNV gain in G3"
    elif delta < -0.25:     verdict = "CNV loss in G3"
    elif delta >  0.10:     verdict = "subtle gain in G3"
    else:                   verdict = "subtle loss in G3"

    print(f"{chrom:10s}  {wt_u:+7.3f}  {wt_a:+7.3f}  {wt_mean:+8.3f}  "
          f"{g3_u:+7.3f}  {g3_a:+7.3f}  {g3_mean:+8.3f}  {delta:+9.3f}  {verdict}")

    redistrib_rows.append({
        "chromosome": chrom,
        "WT-U_log2": round(wt_u, 4), "WT-A_log2": round(wt_a, 4),
        "WT_mean_log2": round(wt_mean, 4),
        "G3-U_log2": round(g3_u, 4), "G3-A_log2": round(g3_a, 4),
        "G3_mean_log2": round(g3_mean, 4),
        "delta_log2": round(delta, 4), "verdict": verdict,
    })

redistrib_df = pd.DataFrame(redistrib_rows)
redistrib_df.to_csv(RESULT_DIR / "cnv_redistribution_chr3_15_18.csv", index=False)
print(f"\nSaved → results/u2os/cnv_redistribution_chr3_15_18.csv\n")

print("\nDone.")
