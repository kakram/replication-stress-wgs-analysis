#!/usr/bin/env python3
"""
Table 5.3.1 generator — MCF-7 somatic variant burden as a Word table.

Reads outputs/ch5_variant_burden.csv and writes two .docx files, each
containing a native Word table that can be copied into the main thesis
document:

  outputs/Table_5_3_1_full.docx       (25 rows — exhaustive reference)
  outputs/Table_5_3_1_condensed.docx  (13 rows — main-text version,
                                       restricted to the metrics named
                                       in the §5.3.1 prose)

Both share the same caption from the §5.3.1 deliverable. Numeric cells
are right-aligned and integer-formatted with thousands separators;
TMB rows render with two decimal places.

Usage:
    python3 scripts/make_ch5_variant_burden_table.py

Requires:
    pip install python-docx pandas
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

try:
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Pt
except ImportError:
    sys.exit(
        "python-docx not installed. Install with:\n"
        "    pip install python-docx"
    )


REPO = Path(__file__).resolve().parent.parent
CSV_IN = REPO / "outputs" / "ch5_variant_burden.csv"
DOCX_FULL = REPO / "outputs" / "Table_5_3_1_full.docx"
DOCX_CONDENSED = REPO / "outputs" / "Table_5_3_1_condensed.docx"

SAMPLE_COLS = ["WTUN", "WTAPH", "B4UN", "B4APH"]
SAMPLE_HEADERS = {
    "WTUN":  "WT untreated",
    "WTAPH": "WT + aphidicolin",
    "B4UN":  "B4 ZFP36L1\u207B\u002F\u207B untreated",
    "B4APH": "B4 ZFP36L1\u207B\u002F\u207B + aphidicolin",
}

CAPTION = (
    "Table 5.3.1. Genome-wide somatic variant burden across four MCF-7 "
    "conditions. Total variants, PASS-filtered counts, filter-class "
    "breakdown, coding consequence distribution, SIFT predictions, "
    "gnomAD novelty, and tumour mutational burden (TMB) are shown for "
    "each sample. Filter labels follow GATK / Mutect2 conventions; rows "
    "carrying multiple filter terms increment each matching filter "
    "column once, so column sums may exceed (Total \u2212 PASS)."
)

# Each entry is (display_label, csv_column, is_subrow_indented)
FULL_ROWS = [
    ("Total variants (VCF)",        "total_variants",            False),
    ("PASS variants",               "pass_variants",             False),
    ("FILTER: germline",            "filter_germline",           True),
    ("FILTER: panel_of_normals",    "filter_panel_of_normals",   True),
    ("FILTER: weak_evidence",       "filter_weak_evidence",      True),
    ("FILTER: clustered_events",    "filter_clustered",          True),
    ("FILTER: other",               "filter_other",              True),
    ("Coding variants (total)",     "coding_total",              False),
    ("Missense",                    "missense",                  True),
    ("Synonymous",                  "synonymous",                True),
    ("Nonsense",                    "nonsense",                  True),
    ("Frameshift deletion",         "frameshift_del",            True),
    ("Frameshift insertion",        "frameshift_ins",            True),
    ("In-frame deletion",           "in_frame_del",              True),
    ("In-frame insertion",          "in_frame_ins",              True),
    ("Splice site",                 "splice_site",               True),
    ("Splice region",               "splice_region",             True),
    ("SIFT deleterious",            "sift_deleterious",          False),
    ("SIFT tolerated",              "sift_tolerated",            False),
    ("SIFT unknown",                "sift_unknown",              False),
    ("gnomAD known (AF > 0.001)",   "gnomad_known",              False),
    ("gnomAD novel (AF \u2264 0.001)", "gnomad_novel",           False),
    ("TMB (mut/Mb)",                "tmb_mut_per_mb",            False),
    ("Genes with \u22651 variant",  "genes_with_any_variant",    False),
    ("Genes with coding variant",   "genes_with_coding_variant", False),
]

# Condensed version — the rows actually referenced in the §5.3.1 prose.
CONDENSED_KEEP = {
    "total_variants",
    "pass_variants",
    "filter_germline",
    "filter_panel_of_normals",
    "filter_weak_evidence",
    "filter_clustered",
    "filter_other",
    "coding_total",
    "missense",
    "sift_deleterious",
    "sift_tolerated",
    "gnomad_novel",
    "tmb_mut_per_mb",
}


def fmt(value, col: str) -> str:
    """Format a cell value: TMB to 2 d.p., everything else as a thousands-
    separated integer."""
    if col == "tmb_mut_per_mb":
        return f"{value:.2f}"
    return f"{int(value):,}"


def build_table(doc, rows, df):
    """Build a Word table from the (label, csv_col, indent) row spec list."""
    n_cols = 1 + len(SAMPLE_COLS)
    table = doc.add_table(rows=1 + len(rows), cols=n_cols)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Light Grid Accent 1"

    # Header row
    hdr = table.rows[0].cells
    hdr[0].text = "Metric"
    for i, s in enumerate(SAMPLE_COLS):
        hdr[i + 1].text = SAMPLE_HEADERS[s]
    for cell in hdr:
        for p in cell.paragraphs:
            for r in p.runs:
                r.bold = True

    # Data rows
    for i, (label, col, indent) in enumerate(rows, start=1):
        row = table.rows[i].cells
        prefix = "    " if indent else ""
        row[0].text = prefix + label
        for j, s in enumerate(SAMPLE_COLS):
            value = df.loc[df["sample"] == s, col].iloc[0]
            row[j + 1].text = fmt(value, col)
            for p in row[j + 1].paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.RIGHT

    # Uniform 10pt font across the table
    for r in table.rows:
        for cell in r.cells:
            for p in cell.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(10)


def write_docx(rows, out_path: Path, df: pd.DataFrame):
    doc = Document()
    cap = doc.add_paragraph()
    cap_run = cap.add_run(CAPTION)
    cap_run.font.size = Pt(10)
    doc.add_paragraph()  # spacer
    build_table(doc, rows, df)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(out_path)
    print(f"  wrote {out_path.relative_to(REPO)}")


def main():
    if not CSV_IN.exists():
        sys.exit(f"Input CSV not found: {CSV_IN}")

    df = pd.read_csv(CSV_IN)

    print("=== Building Table 5.3.1 (full and condensed Word versions) ===\n")
    print(f"Loaded {len(df)} samples from {CSV_IN.relative_to(REPO)}\n")

    write_docx(FULL_ROWS, DOCX_FULL, df)
    condensed = [r for r in FULL_ROWS if r[1] in CONDENSED_KEEP]
    write_docx(condensed, DOCX_CONDENSED, df)

    print("\nDone. Open either .docx in Word, copy the table, paste into the thesis.")


if __name__ == "__main__":
    main()
