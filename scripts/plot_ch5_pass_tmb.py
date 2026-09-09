#!/usr/bin/env python3
"""
Figure 5.3.1 — PASS-variant counts and tumour mutational burden across
the four MCF-7 conditions.

A two-panel bar chart that visualises the genome-wide variant burden
discussed in §5.3.1. Panel A shows PASS-filtered variant counts; Panel B
shows tumour mutational burden (TMB) in mutations per megabase. The
WT-vs-B4 ratio under matched aphidicolin status is annotated above each
genotype pair, making the ~7- to 8-fold asymmetry — and its invariance
under aphidicolin — immediately legible.

Inputs:  outputs/ch5_variant_burden.csv
Outputs:
  figures/mcf7/pass_tmb_burden_mcf7.png
  figures/mcf7/pass_tmb_burden_mcf7.pdf

Usage:
    python3 scripts/plot_ch5_pass_tmb.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


REPO = Path(__file__).resolve().parent.parent
CSV_IN = REPO / "outputs" / "ch5_variant_burden.csv"
FIG_DIR = REPO / "figures" / "mcf7"

SAMPLES = ["WTUN", "WTAPH", "B4UN", "B4APH"]

LABELS = {
    "WTUN":  "WT\nuntreated",
    "WTAPH": "WT\n+ APH",
    "B4UN":  "B4 ZFP36L1\u207B\u002F\u207B\nuntreated",
    "B4APH": "B4 ZFP36L1\u207B\u002F\u207B\n+ APH",
}

# Cool palette for WT (parental), warm palette for B4 (knockout).
# Treatment darkens the same hue.
COLOURS = {
    "WTUN":  "#9EC1DC",
    "WTAPH": "#2C5F8D",
    "B4UN":  "#F2B36A",
    "B4APH": "#B85A1A",
}


def annotate_ratio(ax, x1: int, y1: float, x2: int, y2: float, ratio_text: str):
    """Draw a thin bracket between two bars and label it with a ratio."""
    y_max = max(y1, y2)
    pad = y_max * 0.08
    bracket_y = y_max + pad
    ax.plot([x1, x1, x2, x2],
            [y_max + pad * 0.4, bracket_y, bracket_y, y_max + pad * 0.4],
            color="#333333", linewidth=1)
    ax.text((x1 + x2) / 2, bracket_y + pad * 0.2, ratio_text,
            ha="center", va="bottom", fontsize=9, color="#333333")


def panel_bar(ax, values, *, title: str, ylabel: str, value_fmt: str,
              ratio_fmt: str):
    """Render one of the two panels with consistent styling."""
    x = list(range(len(SAMPLES)))
    colours = [COLOURS[s] for s in SAMPLES]
    labels = [LABELS[s] for s in SAMPLES]

    ax.bar(x, values, color=colours, edgecolor="#222", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_ylim(0, max(values) * 1.28)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)

    for i, v in enumerate(values):
        ax.text(i, v * 1.01, value_fmt.format(v), ha="center", va="bottom",
                fontsize=8, color="#222")

    r_un = values[0] / values[2]
    r_aph = values[1] / values[3]
    annotate_ratio(ax, 0, values[0], 2, values[2], ratio_fmt.format(r_un))
    annotate_ratio(ax, 1, values[1], 3, values[3], ratio_fmt.format(r_aph))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    if not CSV_IN.exists():
        sys.exit(f"Input CSV not found: {CSV_IN}")

    df = pd.read_csv(CSV_IN).set_index("sample").loc[SAMPLES]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 4.6))

    panel_bar(
        ax_a, df["pass_variants"].values,
        title="A. PASS-filtered variant count",
        ylabel="PASS variants (count)",
        value_fmt="{:,.0f}",
        ratio_fmt="ratio = {:.2f}\u00D7",
    )
    panel_bar(
        ax_b, df["tmb_mut_per_mb"].values,
        title="B. Tumour mutational burden",
        ylabel="TMB (mutations / Mb)",
        value_fmt="{:.2f}",
        ratio_fmt="ratio = {:.2f}\u00D7",
    )

    fig.tight_layout()

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / "pass_tmb_burden_mcf7.png"
    pdf_path = FIG_DIR / "pass_tmb_burden_mcf7.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path.relative_to(REPO)}")
    print(f"Wrote {pdf_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
