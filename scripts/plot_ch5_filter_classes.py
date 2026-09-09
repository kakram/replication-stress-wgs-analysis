#!/usr/bin/env python3
"""
Figure 5.3.2 — Filter-class composition across the four MCF-7 conditions.

Each TNfilter rejection class is expressed as a percentage of total variants
called per sample. Four bars per class (one per sample) make the
within-genotype conservation (WT-UN vs WT-APH; B4-UN vs B4-APH) and the
cross-genotype divergence (WT vs B4) directly comparable. The
panel_of_normals class is zero in all four samples and is omitted from
the plot — that fact is noted in the §5.3.1 prose and the figure caption.

Filter classes can overlap (a single variant may be flagged by more than
one filter), so the bar heights for each sample do not sum to 100 %.
The y-axis is the prevalence of each filter class within the per-sample
variant call set, not a partition.

Inputs:  outputs/ch5_variant_burden.csv
Outputs:
  figures/mcf7/filter_class_proportions_mcf7.png
  figures/mcf7/filter_class_proportions_mcf7.pdf

Usage:
    python3 scripts/plot_ch5_filter_classes.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parent.parent
CSV_IN = REPO / "outputs" / "ch5_variant_burden.csv"
FIG_DIR = REPO / "figures" / "mcf7"

SAMPLES = ["WTUN", "WTAPH", "B4UN", "B4APH"]
SAMPLE_LABELS = {
    "WTUN":  "WT untreated",
    "WTAPH": "WT + APH",
    "B4UN":  "B4 untreated",
    "B4APH": "B4 + APH",
}
SAMPLE_COLOURS = {
    "WTUN":  "#9EC1DC",
    "WTAPH": "#2C5F8D",
    "B4UN":  "#F2B36A",
    "B4APH": "#B85A1A",
}

# (display label, csv column). panel_of_normals omitted: zero in all samples.
CLASSES = [
    ("PASS",             "pass_variants"),
    ("germline",         "filter_germline"),
    ("weak_evidence",    "filter_weak_evidence"),
    ("clustered_events", "filter_clustered"),
    ("other",            "filter_other"),
]


def main():
    if not CSV_IN.exists():
        sys.exit(f"Input CSV not found: {CSV_IN}")

    df = pd.read_csv(CSV_IN).set_index("sample").loc[SAMPLES]

    # Build a (n_classes x n_samples) matrix of percentages of total variants.
    pct = np.zeros((len(CLASSES), len(SAMPLES)))
    for j, s in enumerate(SAMPLES):
        total = df.loc[s, "total_variants"]
        for i, (_, col) in enumerate(CLASSES):
            pct[i, j] = 100.0 * df.loc[s, col] / total

    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    n_classes, n_samples = pct.shape
    bar_w = 0.18
    x_base = np.arange(n_classes)

    for j, s in enumerate(SAMPLES):
        offset = (j - (n_samples - 1) / 2) * bar_w
        bars = ax.bar(x_base + offset, pct[:, j], width=bar_w,
                      color=SAMPLE_COLOURS[s], edgecolor="#222",
                      linewidth=0.5, label=SAMPLE_LABELS[s])
        for k, b in enumerate(bars):
            ax.text(b.get_x() + b.get_width() / 2,
                    b.get_height() + max(pct.flatten()) * 0.012,
                    f"{pct[k, j]:.1f}",
                    ha="center", va="bottom",
                    fontsize=7, color="#222")

    ax.set_xticks(x_base)
    ax.set_xticklabels([label for label, _ in CLASSES], fontsize=10)
    ax.set_ylabel("Percentage of total variants (%)", fontsize=10)
    ax.set_title(
        "Filter-class prevalence per sample "
        "(percentage of total variants; classes may overlap)",
        fontsize=11, loc="left",
    )
    ax.set_ylim(0, max(pct.flatten()) * 1.18)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper right", ncol=2)

    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / "filter_class_proportions_mcf7.png"
    pdf_path = FIG_DIR / "filter_class_proportions_mcf7.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path.relative_to(REPO)}")
    print(f"Wrote {pdf_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
