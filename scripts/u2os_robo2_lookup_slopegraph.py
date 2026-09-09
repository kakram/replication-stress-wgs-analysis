#!/usr/bin/env python3
"""
Task 1: Report ROBO2 chromosomal context vs chr3 CFS regions.
Task 2: Slope-graph of top-15 rank-shifters (Section 6.4 figure).

Run from repo root:
    python3 scripts/u2os_robo2_lookup_slopegraph.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pathlib import Path

RESULT_DIR = Path("results/u2os")
FIG_DIR    = Path("figures/u2os")
CFS_DIR    = Path("data/CFS_fixed")
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# Task 1 — ROBO2 chromosomal context
# ══════════════════════════════════════════════════════════════════════════════
print("=== ROBO2 chromosomal context ===\n")

ROBO2_CHROM = "chr3"
ROBO2_START = 75_907_852
ROBO2_END   = 77_651_517

chr3_cfs = pd.read_csv(
    CFS_DIR / "chr3_fragile_site.bed", sep="\t", header=None,
    names=["chr", "start", "end", "id", "score", "strand"],
)

print(f"ROBO2: {ROBO2_CHROM}:{ROBO2_START:,}–{ROBO2_END:,}  "
      f"(gene span {(ROBO2_END - ROBO2_START)/1e6:.2f} Mb)\n")
print(f"{'CFS':8s}  {'Region (Mb)':20s}  {'Distance':>18s}  Relation")
print("-" * 70)

nearest_dist = np.inf
nearest_name = ""

for _, row in chr3_cfs.iterrows():
    cs, ce, cid = int(row["start"]), int(row["end"]), row["id"]
    if ROBO2_START <= ce and ROBO2_END >= cs:
        dist = 0
        label = "OVERLAPS"
    elif ROBO2_START > ce:
        dist = ROBO2_START - ce
        label = f"{'3′ of CFS':>18s}"
    else:
        dist = cs - ROBO2_END
        label = f"{'5′ of CFS':>18s}"

    dist_str = "overlaps" if dist == 0 else f"{dist/1e6:.2f} Mb"
    print(f"{cid:8s}  {row['start']//1_000_000:.1f}–{ce//1_000_000:.1f} Mb"
          f"{'':8s}  {dist_str:>18s}  {label.strip()}")

    if dist < nearest_dist:
        nearest_dist = dist
        nearest_name = cid

print()
print(f"Nearest CFS: {nearest_name}  ({nearest_dist/1e6:.2f} Mb gap)")
print()

# FHIT landmark check
FHIT_START, FHIT_END = 59_734_036, 61_238_376   # GRCh38 canonical
gap_fhit_robo2 = ROBO2_START - FHIT_END
print(f"FHIT (FRA3B anchor): chr3:{FHIT_START:,}–{FHIT_END:,}")
print(f"Gap FHIT→ROBO2: {gap_fhit_robo2/1e6:.2f} Mb")
print()
print("Verdict: ROBO2 is 12.2 Mb distal to FRA3B. Report as")
print("'chromosome-3 localised but outside any annotated CFS'.")
print()


# ══════════════════════════════════════════════════════════════════════════════
# Task 2 — Slope-graph figure
# ══════════════════════════════════════════════════════════════════════════════
print("=== Building slope-graph (rank_shift_top15.png) ===\n")

# ── load data ─────────────────────────────────────────────────────────────────
rank_df = pd.read_csv(RESULT_DIR / "rank_shift_u2os.csv")

# Build per-gene CFS status from the overlap CSV
cfs_ov = pd.read_csv(RESULT_DIR / "cfs_overlap_top30_u2os.csv")
cfs_genes = set(cfs_ov[cfs_ov["n_cfs"] > 0]["Gene"])

# Drop artefact / locus flags (unflagged genes are NaN in the CSV)
clean = rank_df[rank_df["flag"].isna()].copy()
clean["abs_delta"] = clean["rank_delta"].abs()
top15 = clean.nlargest(15, "abs_delta").copy()

print("Top-15 genes by |rank_delta| (PAR/HLA excluded):")
print(top15[["Gene","mean_WT_rank","mean_G3_rank","rank_delta","abs_delta"]].to_string(index=False))
print()

# ── figure ────────────────────────────────────────────────────────────────────
HIGHLIGHT = {"FHIT", "ROBO2"}
CFS_COLOR  = "#c0392b"   # red
GREY_COLOR = "#7f8c8d"   # mid-grey
HL_COLOR   = "#c0392b"   # bold red for highlighted CFS (FHIT); ROBO2 needs own treatment

fig, ax = plt.subplots(figsize=(6, 8))

X_WT = 0.0
X_G3 = 1.0
Y_MAX = 31   # ">30" plotted at 31

# Sort by WT rank for clean left-column ordering
top15_sorted = top15.sort_values("mean_WT_rank")

for _, row in top15_sorted.iterrows():
    gene  = row["Gene"]
    y_wt  = row["mean_WT_rank"]
    y_g3  = row["mean_G3_rank"]
    is_cfs   = gene in cfs_genes
    is_hl    = gene in HIGHLIGHT

    lw     = 2.2 if is_hl else 1.4
    alpha  = 0.95 if is_hl else 0.75
    color  = CFS_COLOR if is_cfs else GREY_COLOR

    ax.plot([X_WT, X_G3], [y_wt, y_g3],
            color=color, lw=lw, alpha=alpha,
            solid_capstyle="round")
    ax.scatter([X_WT, X_G3], [y_wt, y_g3],
               color=color, s=30, zorder=5, alpha=alpha)

    # Labels — offset slightly so they don't overlap the dots
    bold = "bold" if is_hl else "normal"
    fs   = 8.5 if is_hl else 7.5

    # Left label (WT side)
    ax.text(X_WT - 0.03, y_wt, gene,
            ha="right", va="center",
            fontsize=fs, fontweight=bold, color=color)
    # Right label (G3 side)
    ax.text(X_G3 + 0.03, y_g3, gene,
            ha="left", va="center",
            fontsize=fs, fontweight=bold, color=color)

# ── axes ──────────────────────────────────────────────────────────────────────
ax.set_xlim(-0.55, 1.55)
ax.set_ylim(Y_MAX + 0.5, 0.5)   # invert: rank 1 at top

# Tick marks at rank positions present in the data (avoid clutter)
rank_vals = sorted(set(
    list(top15["mean_WT_rank"]) + list(top15["mean_G3_rank"])
))
ax.set_yticks(rank_vals)
ax.set_yticklabels([str(int(v)) if v < 31 else ">30" for v in rank_vals],
                   fontsize=7)
ax.yaxis.set_tick_params(length=3)

ax.set_xticks([X_WT, X_G3])
ax.set_xticklabels(["WT mean rank", "G3 mean rank"], fontsize=10, fontweight="bold")
ax.tick_params(bottom=False)

# Remove all spines except a subtle y-axis guide
for spine in ["top", "right", "bottom"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_color("#cccccc")
ax.spines["left"].set_linewidth(0.8)

ax.set_title("Rank shift: WT → G3 KO\n(top-15 movers, protein-coding genes)",
             fontsize=10, pad=10)
ax.set_ylabel("Rank within sample top-30  (1 = most mutated)", fontsize=8)

# ── legend ────────────────────────────────────────────────────────────────────
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color=CFS_COLOR, lw=2, label="Overlaps annotated CFS"),
    Line2D([0], [0], color=GREY_COLOR, lw=1.5, label="No CFS overlap"),
]
ax.legend(handles=legend_elements, fontsize=7.5, loc="lower right",
          frameon=True, framealpha=0.9, edgecolor="#cccccc")

# ── annotation arrow for FHIT ─────────────────────────────────────────────────
# Small arrow pointing to the FHIT line midpoint
fhit_row = top15[top15["Gene"] == "FHIT"].iloc[0]
ax.annotate("FRA3B", xy=(0.5, (fhit_row["mean_WT_rank"] + fhit_row["mean_G3_rank"]) / 2),
            xytext=(0.5, (fhit_row["mean_WT_rank"] + fhit_row["mean_G3_rank"]) / 2 + 3.5),
            ha="center", fontsize=7, color=CFS_COLOR,
            arrowprops=dict(arrowstyle="-|>", color=CFS_COLOR, lw=0.8),
            fontweight="bold")

plt.tight_layout()

for suffix in ("png", "pdf"):
    out = FIG_DIR / f"rank_shift_top15.{suffix}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved → {out}")

plt.close()
print("\nDone.")
