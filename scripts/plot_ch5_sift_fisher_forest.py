#!/usr/bin/env python3
"""
Figure 5.10.1 — Forest plot of SIFT del:tol Fisher exact test odds ratios.

Visualises the three pairwise Fisher's-exact test results from §5.10 as
a forest plot on a log-scaled odds-ratio axis. Each test contributes one
row showing the point estimate (filled circle for non-significant
within-genotype contrasts; filled diamond for the between-genotype
contrast that survives Bonferroni correction) with a horizontal 95 %
confidence interval bar capped at each end. A vertical dashed reference
line at OR = 1 marks the null hypothesis.

A right-hand text panel reports the OR, 95 % CI, and p-value for each
row, with the Bonferroni-significant result highlighted in bold colour.

Inputs:  outputs/ch5_sift_fisher_results.csv
Outputs:
  figures/mcf7/sift_fisher_forest_mcf7.png
  figures/mcf7/sift_fisher_forest_mcf7.pdf

Usage:
    python3 scripts/plot_ch5_sift_fisher_forest.py
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
CSV_IN = REPO / "outputs" / "ch5_sift_fisher_results.csv"
FIG_DIR = REPO / "figures" / "mcf7"

# Top-to-bottom order of rows in the plot (visually first row at the top)
ROW_ORDER = [
    "Test 1: WTUN vs WTAPH (within-WT treatment)",
    "Test 2: B4UN vs B4APH (within-B4 treatment)",
    "Test 3: WT pooled vs B4 pooled (between-genotype)",
]

ROW_LABELS = {
    "Test 1: WTUN vs WTAPH (within-WT treatment)":
        "WT untreated\nvs WT + APH",
    "Test 2: B4UN vs B4APH (within-B4 treatment)":
        "B4 untreated\nvs B4 + APH",
    "Test 3: WT pooled vs B4 pooled (between-genotype)":
        "WT pooled\nvs B4 pooled",
}

ALPHA_CORR = 0.05 / 3

COLOUR_NONSIG = "#7C8FA1"  # muted slate-blue
COLOUR_SIG    = "#1F4E79"  # deep navy


def format_p(p: float) -> str:
    """Format p-value: scientific notation if very small, decimal otherwise."""
    if p < 0.001:
        exponent = int(np.floor(np.log10(p)))
        mantissa = p / (10 ** exponent)
        super_map = str.maketrans(
            "-0123456789",
            "\u207B\u2070\u00B9\u00B2\u00B3\u2074\u2075\u2076\u2077\u2078\u2079",
        )
        exp_str = str(exponent).translate(super_map)
        return f"{mantissa:.1f} \u00D7 10{exp_str}"
    return f"{p:.2f}"


def main():
    if not CSV_IN.exists():
        sys.exit(
            f"Fisher results CSV not found: {CSV_IN}\n"
            f"Run scripts/ch5_sift_fisher_test.py first to generate it."
        )

    df = pd.read_csv(CSV_IN).set_index("test").loc[ROW_ORDER].reset_index()
    n_tests = len(df)

    fig, (ax_plot, ax_text) = plt.subplots(
        1, 2, figsize=(11.5, 4.0),
        gridspec_kw={"width_ratios": [2.4, 1.6]},
    )

    # Y-positions: top-to-bottom (first row at the top of the plot)
    y_positions = list(range(n_tests))[::-1]

    # ── Forest plot panel ────────────────────────────────────────────────
    for i, fr in enumerate(df.itertuples(index=False)):
        y = y_positions[i]
        is_sig = fr.p_value < ALPHA_CORR
        colour = COLOUR_SIG if is_sig else COLOUR_NONSIG
        marker = "D" if is_sig else "o"
        marker_size = 110 if is_sig else 75

        # CI horizontal bar
        ax_plot.plot(
            [fr.ci95_low, fr.ci95_high], [y, y],
            color=colour, linewidth=2.0, solid_capstyle="round", zorder=2,
        )
        # End caps on the CI bar
        for x_end in (fr.ci95_low, fr.ci95_high):
            ax_plot.plot(
                [x_end, x_end], [y - 0.14, y + 0.14],
                color=colour, linewidth=2.0, zorder=2,
            )
        # Point estimate marker
        ax_plot.scatter(
            [fr.odds_ratio], [y], s=marker_size, color=colour,
            marker=marker, zorder=3, edgecolor="white", linewidth=0.8,
        )

    # Vertical reference line at OR = 1
    ax_plot.axvline(1.0, linestyle="--", color="#888888",
                    linewidth=1.0, zorder=1)

    # Y-axis labels
    ax_plot.set_yticks(y_positions)
    ax_plot.set_yticklabels(
        [ROW_LABELS[t] for t in df["test"]],
        fontsize=10,
    )

    # X-axis: log scale, range chosen to contain all CIs with headroom
    ax_plot.set_xscale("log")
    ax_plot.set_xlim(0.25, 3.0)
    ax_plot.set_xticks([0.3, 0.5, 1.0, 2.0, 3.0])
    ax_plot.set_xticklabels(["0.3", "0.5", "1.0", "2.0", "3.0"])
    ax_plot.set_xlabel("Odds ratio  (Group 1 vs Group 2, log scale)",
                       fontsize=10)
    ax_plot.set_title(
        "SIFT del:tol Fisher's-exact tests\u2014odds ratios with 95 % CIs",
        fontsize=11, loc="left", pad=10,
    )

    # Cosmetic
    ax_plot.spines["top"].set_visible(False)
    ax_plot.spines["right"].set_visible(False)
    ax_plot.set_ylim(-0.6, n_tests - 0.4)
    ax_plot.grid(axis="x", which="major", linestyle=":", alpha=0.4)
    ax_plot.set_axisbelow(True)

    # ── Text annotation panel ────────────────────────────────────────────
    ax_text.set_xlim(0, 1)
    ax_text.set_ylim(-0.6, n_tests - 0.4)
    ax_text.axis("off")

    for i, fr in enumerate(df.itertuples(index=False)):
        y = y_positions[i]
        is_sig = fr.p_value < ALPHA_CORR
        colour = COLOUR_SIG if is_sig else COLOUR_NONSIG
        weight = "bold" if is_sig else "normal"

        or_line = (f"OR = {fr.odds_ratio:.2f}  "
                   f"(95 % CI: {fr.ci95_low:.2f}\u2013{fr.ci95_high:.2f})")
        p_line = f"p = {format_p(fr.p_value)}"
        if is_sig:
            p_line += "  \u2605 significant"
        else:
            p_line += "    n.s."

        annotation = f"{or_line}\n{p_line}"
        ax_text.text(
            0.02, y, annotation,
            va="center", ha="left",
            fontsize=9, color=colour, weight=weight,
            family="monospace",
        )

    fig.tight_layout()

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / "sift_fisher_forest_mcf7.png"
    pdf_path = FIG_DIR / "sift_fisher_forest_mcf7.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {png_path.relative_to(REPO)}")
    print(f"Wrote {pdf_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
