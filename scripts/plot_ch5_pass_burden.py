#!/usr/bin/env python3
"""
Figure 5.X (formerly 5.3.1) — PASS-variant counts across the four MCF-7
conditions. Single panel.

Supersedes plot_ch5_pass_tmb.py. The former panel B (provider-reported
tumour mutational burden) was dropped on 2026-09-09: the provider's TMB
values (0.05–0.42 /Mb) cannot be reproduced from the PASS counts over the
978.35 Mb callable genome (~147 /Mb), so the variant subset and
denominator behind them are unknown and the figure cannot be defended.

The panel shows PASS-filtered small-variant counts for WT and B4, each
untreated and aphidicolin-treated, with the WT-to-B4 ratio under matched
treatment annotated above each genotype pair. The calling design that
produced each call set (tumour-only for WT, paired against the
treatment-matched WT for B4) is labelled beneath the bars, because that
design — not biology — is the principal cause of the ~7- to 8-fold gap.

Inputs:  outputs/ch5_variant_burden.csv
Outputs:
  figures/mcf7/pass_burden_mcf7.png
  figures/mcf7/pass_burden_mcf7.pdf

Usage:
    python3 scripts/plot_ch5_pass_burden.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea, VPacker
import pandas as pd


REPO = Path(__file__).resolve().parent.parent
CSV_IN = REPO / "outputs" / "ch5_variant_burden.csv"
FIG_DIR = REPO / "figures" / "mcf7"

SAMPLES = ["WTUN", "WTAPH", "B4UN", "B4APH"]

# Tick labels are built from text pieces so the gene symbol can be
# italic while "B4" and the genotype superscript stay upright. mathtext
# cannot do this: it forces digits to roman inside \mathit.
LABELS = {
    "WTUN":  [("WT", {})],
    "WTAPH": [("WT", {})],
    "B4UN":  [("B4 ", {}), ("ZFP36L1", {"style": "italic"}),
              ("\u207B\u141F\u207B", {})],
    "B4APH": [("B4 ", {}), ("ZFP36L1", {"style": "italic"}),
              ("\u207B\u141F\u207B", {})],
}
TREATMENT = {
    "WTUN": "untreated", "WTAPH": "+ APH",
    "B4UN": "untreated", "B4APH": "+ APH",
}

# Cool palette for WT (parental), warm palette for B4 (knockout).
# Treatment darkens the same hue.
COLOURS = {
    "WTUN":  "#9EC1DC",
    "WTAPH": "#2C5F8D",
    "B4UN":  "#F2B36A",
    "B4APH": "#B85A1A",
}

# Calling design per genotype group, drawn beneath the x-axis
DESIGN = [
    (0, 1, "tumour-only call"),
    (2, 3, "paired call vs treatment-matched WT"),
]


def annotate_ratio(ax, x1: int, x2: int, bracket_y: float, tick: float,
                   text: str):
    """Thin bracket between two bars at an explicit height, labelled with
    the ratio. Heights are passed in so the two brackets can be staggered
    clear of each other and of the bar-value labels."""
    ax.plot([x1, x1, x2, x2],
            [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
            color="#333333", linewidth=1)
    ax.text((x1 + x2) / 2, bracket_y + tick * 0.3, text,
            ha="center", va="bottom", fontsize=9, color="#333333")


def tick_label(ax, x: int, pieces, second_line: str, fontsize: float = 9):
    """Two-line tick label: a first line composed of styled pieces, a
    plain second line. Anchored just below the axis at bar position x."""
    line1 = HPacker(children=[
        TextArea(text, textprops=dict(fontsize=fontsize, **props))
        for text, props in pieces], align="baseline", pad=0, sep=0)
    line2 = TextArea(second_line, textprops=dict(fontsize=fontsize))
    box = VPacker(children=[line1, line2], align="center", pad=0, sep=2)
    ab = AnnotationBbox(box, (x, 0), xybox=(0, -6), xycoords=("data", "axes fraction"),
                        boxcoords="offset points", box_alignment=(0.5, 1.0),
                        frameon=False, pad=0)
    ax.add_artist(ab)


def annotate_design(ax, x1: int, x2: int, text: str, y: float):
    """Bracket beneath a group of bars, labelled with the calling design."""
    ax.plot([x1 - 0.35, x1 - 0.35, x2 + 0.35, x2 + 0.35],
            [y + 0.012, y, y, y + 0.012],
            color="#555555", linewidth=0.8, clip_on=False,
            transform=ax.get_xaxis_transform())
    ax.text((x1 + x2) / 2, y - 0.02, text, ha="center", va="top",
            fontsize=8, color="#555555", style="italic",
            transform=ax.get_xaxis_transform())


def main():
    if not CSV_IN.exists():
        sys.exit(f"Input CSV not found: {CSV_IN}")

    df = pd.read_csv(CSV_IN).set_index("sample").loc[SAMPLES]
    values = df["pass_variants"].values

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    x = list(range(len(SAMPLES)))

    ax.bar(x, values, color=[COLOURS[s] for s in SAMPLES],
           edgecolor="#222", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([""] * len(SAMPLES))
    ax.tick_params(axis="x", length=3)
    for i, smp in enumerate(SAMPLES):
        tick_label(ax, i, LABELS[smp], TREATMENT[smp])
    ax.set_ylabel("PASS-filtered variants (count)", fontsize=10)
    top = max(values)
    ax.set_ylim(0, top * 1.42)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))

    for i, v in enumerate(values):
        ax.text(i, v * 1.01, f"{v:,.0f}", ha="center", va="bottom",
                fontsize=8, color="#222")

    tick = top * 0.03
    annotate_ratio(ax, 0, 2, top * 1.12, tick,
                   f"ratio = {values[0] / values[2]:.2f}\u00D7")
    annotate_ratio(ax, 1, 3, top * 1.26, tick,
                   f"ratio = {values[1] / values[3]:.2f}\u00D7")

    for x1, x2, text in DESIGN:
        annotate_design(ax, x1, x2, text, y=-0.17)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.26)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / "pass_burden_mcf7.png"
    pdf_path = FIG_DIR / "pass_burden_mcf7.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path.relative_to(REPO)}")
    print(f"Wrote {pdf_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
