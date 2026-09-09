#!/usr/bin/env python3
"""
Table 5.10.1 generator — SIFT del:tol Fisher exact test results as a Word table.

Reads outputs/ch5_sift_fisher_results.csv (the three pairwise Fisher tests)
and outputs/ch5_variant_burden.csv (used to recompute the 2x4 chi-square
omnibus test for the table footnote), and writes a .docx file containing
a native Word table ready to copy into the §5.10 prose.

Outputs:
  outputs/Table_5_10_1.docx

The table summarises three pairwise Fisher tests (within-WT treatment,
within-B4 treatment, between-genotype pooled across treatment) with
effect sizes, exact 95 % confidence intervals, two-sided p-values, and
Bonferroni-coded significance. The 2x4 chi-square omnibus test is
reported as a table footnote.

Usage:
    python3 scripts/make_ch5_sift_fisher_table.py

Requires:
    pip install python-docx pandas scipy
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

try:
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Pt
except ImportError:
    sys.exit("python-docx not installed. Install with: pip install python-docx")


REPO = Path(__file__).resolve().parent.parent
FISHER_CSV = REPO / "outputs" / "ch5_sift_fisher_results.csv"
BURDEN_CSV = REPO / "outputs" / "ch5_variant_burden.csv"
DOCX_OUT = REPO / "outputs" / "Table_5_10_1.docx"

# Short, clean row labels for the table — the CSV test column carries the
# longer analytical labels used in the markdown summary.
ROW_LABELS = {
    "Test 1: WTUN vs WTAPH (within-WT treatment)":        "WT: untreated vs aphidicolin",
    "Test 2: B4UN vs B4APH (within-B4 treatment)":        "B4: untreated vs aphidicolin",
    "Test 3: WT pooled vs B4 pooled (between-genotype)":  "WT vs B4 (pooled across treatment)",
}

ALPHA_CORR = 0.05 / 3  # Bonferroni-corrected significance threshold

CAPTION = (
    "Table 5.10.1. Fisher's exact tests of the ratio of SIFT-deleterious to "
    "SIFT-tolerated missense calls across the four MCF-7 conditions. Three "
    "pairwise tests are reported, with Bonferroni correction across the family "
    f"(\u03B1_corr = 0.05 / 3 = {ALPHA_CORR:.4f}). Counts shown as "
    "(deleterious : tolerated, % deleterious). Significance codes: *** "
    "p < 0.001, ** p < 0.01, * p < 0.05 (all after Bonferroni correction); "
    "n.s. = not significant."
)


def format_p(p: float) -> str:
    """Format a p-value: scientific notation if very small, decimal otherwise."""
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


def significance_code(p: float, alpha_corr: float) -> str:
    """Return a thesis-style significance string for a p-value."""
    if p >= alpha_corr:
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    return "*"


def compute_omnibus_chi2() -> tuple[float, int, float]:
    """Recompute the 2x4 chi-square (sample x SIFT class) for the footnote."""
    df = pd.read_csv(BURDEN_CSV).set_index("sample")
    samples = ("WTUN", "WTAPH", "B4UN", "B4APH")
    table = np.array([
        [df.loc[s, "sift_deleterious"] for s in samples],
        [df.loc[s, "sift_tolerated"] for s in samples],
    ])
    chi2, p, dof, _ = chi2_contingency(table)
    return float(chi2), int(dof), float(p)


def build_table(doc, fisher_df: pd.DataFrame) -> None:
    """Build the Fisher test results table in the document."""
    headers = [
        "Contrast",
        "Group 1 (del : tol, % del)",
        "Group 2 (del : tol, % del)",
        "OR (95 % CI)",
        "p-value",
        "Sig.",
    ]
    table = doc.add_table(rows=1 + len(fisher_df), cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Light Grid Accent 1"

    # Header row — bold
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = h
        for p in cell.paragraphs:
            for r in p.runs:
                r.bold = True
                r.font.size = Pt(10)

    # Data rows
    for row_idx, fr in enumerate(fisher_df.itertuples(index=False), start=1):
        label = ROW_LABELS.get(fr.test, fr.test)
        g1 = f"{fr.g1_del} : {fr.g1_tol} ({fr.g1_pct_del:.1f} %)"
        g2 = f"{fr.g2_del} : {fr.g2_tol} ({fr.g2_pct_del:.1f} %)"
        or_text = (f"{fr.odds_ratio:.2f} "
                   f"({fr.ci95_low:.2f}\u2013{fr.ci95_high:.2f})")
        p_text = format_p(fr.p_value)
        sig = significance_code(fr.p_value, ALPHA_CORR)

        cells = table.rows[row_idx].cells
        cells[0].text = label
        cells[1].text = g1
        cells[2].text = g2
        cells[3].text = or_text
        cells[4].text = p_text
        cells[5].text = sig

        # Bold the entire row if the result survives Bonferroni correction
        is_significant = fr.p_value < ALPHA_CORR
        for c in cells:
            for p in c.paragraphs:
                for r in p.runs:
                    r.bold = is_significant
                    r.font.size = Pt(10)

        # Centre-align numeric / counts columns, left-align label
        for i, c in enumerate(cells):
            for p in c.paragraphs:
                p.alignment = (WD_ALIGN_PARAGRAPH.LEFT if i == 0
                               else WD_ALIGN_PARAGRAPH.CENTER)


def main():
    if not FISHER_CSV.exists():
        sys.exit(
            f"Fisher results CSV not found: {FISHER_CSV}\n"
            f"Run scripts/ch5_sift_fisher_test.py first to generate it."
        )
    if not BURDEN_CSV.exists():
        sys.exit(f"Burden CSV not found: {BURDEN_CSV}")

    fisher_df = pd.read_csv(FISHER_CSV)
    chi2, dof, p_chi2 = compute_omnibus_chi2()

    print("=== Building Table 5.10.1 ===\n")
    print(f"Loaded {len(fisher_df)} pairwise tests "
          f"from {FISHER_CSV.relative_to(REPO)}")
    print(f"Omnibus chi-square: \u03C7\u00B2 = {chi2:.2f}, "
          f"dof = {dof}, p = {p_chi2:.4g}\n")

    doc = Document()

    # Caption
    cap = doc.add_paragraph()
    cap.add_run(CAPTION).font.size = Pt(10)
    doc.add_paragraph()  # spacer

    # The table itself
    build_table(doc, fisher_df)
    doc.add_paragraph()  # spacer

    # Footnote
    fn = doc.add_paragraph()
    fn_text = (
        f"Omnibus 2 \u00D7 4 chi-square (sample \u00D7 SIFT class) confirms global "
        f"heterogeneity: \u03C7\u00B2 = {chi2:.2f}, dof = {dof}, "
        f"p = {format_p(p_chi2)}. Odds ratios and exact 95 % confidence "
        "intervals computed via scipy.stats.contingency.odds_ratio."
    )
    fn_run = fn.add_run(fn_text)
    fn_run.italic = True
    fn_run.font.size = Pt(9)

    DOCX_OUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(DOCX_OUT)
    print(f"Wrote {DOCX_OUT.relative_to(REPO)}")
    print("\nDone. Open Table_5_10_1.docx in Word, copy the table and footnote, "
          "paste into §5.10.")


if __name__ == "__main__":
    main()
