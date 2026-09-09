#!/usr/bin/env python3
"""
Filter-cascade tables for Chapter 5 (Part C of the 9 Sept handover).

Reads the per-sample cascade JSONs written by scripts/qc/filter_cascade.py
(primary-contig set by default) and the BAM depth summary from
scripts/qc/bam_metrics.py, and writes one Word file with three captioned
tables ready to paste into the thesis:

  Table A  Variant counts at each filter tier, with Ts/Tv at PASS and at
           the final tier.
  Table B  Depth-matched treatment contrast: untreated vs treated at the
           final tier for each genotype, alongside sequencing depth and,
           for the paired B4 sets, the tumour/normal depth ratio.
  Table C  Clonal bottleneck diagnostic: alt-read support in the matched
           wild-type among T4 variants of the B4 sets.

Inputs:
  results/cascade_primary/<SAMPLE>_cascade.json  (4 MCF-7 samples)
  results/qc/bam_metrics_summary.tsv
Outputs:
  outputs/Table_5_X_filter_cascade.docx

Usage:
  python3 scripts/make_ch5_filter_cascade_table.py [--cascade-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

try:
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt
except ImportError:
    sys.exit("python-docx not installed: pip install python-docx")

REPO = Path(__file__).resolve().parent.parent
CASCADE_DIR = REPO / "results" / "cascade_primary"
BAM_TSV = REPO / "results" / "qc" / "bam_metrics_summary.tsv"
DOCX_OUT = REPO / "outputs" / "Table_5_X_filter_cascade.docx"

SAMPLES = ["WTUN", "WTAPH", "B4UN", "B4APH"]
LABELS = {
    "WTUN":  "WT untreated",
    "WTAPH": "WT + APH",
    "B4UN":  "B4 untreated",
    "B4APH": "B4 + APH",
}
TIERS = ["T0_ALL", "T1_PASS", "T2_DEPTH", "T3_VAF", "T4_POP", "T5_PRIVATE"]
TIER_LABELS = ["All records", "PASS", "Depth", "VAF", "Population", "Private"]
TIER_HEADERS = ["All", "PASS", "Depth", "VAF", "Pop. AF", "Private"]

CAPTION_A = (
    "Table 5.X. Somatic small-variant counts across a cumulative filter "
    "cascade, MCF-7, primary contigs. Tiers are cumulative: PASS, the "
    "caller's own filter; Depth, read depth \u2265 20 in the tumour and, for "
    "paired call sets, in the matched normal; VAF, allele fraction \u2265 0.05 "
    "with \u2265 3 supporting reads; Population, gnomAD population allele "
    "frequency \u2264 10\u207B\u00B3; Private, no alt-supporting read in the "
    "matched normal (paired call sets only). Wild-type samples were called "
    "in tumour-only mode and have no Private tier. Ts/Tv is the "
    "transition-to-transversion ratio among single-nucleotide variants at "
    "the PASS tier and at the final tier for each sample."
)
CAPTION_B = (
    "Table 5.X. Treatment contrast at the final filter tier, with "
    "sequencing depth. Counts are at the Population tier for wild-type "
    "(tumour-only) and the Private tier for B4 (paired). Mean coverage is "
    "post-deduplication, computed from the BAM files as in Table 5.X. For "
    "the paired B4 call sets the tumour/normal depth ratio governs calling "
    "power in both directions: deeper tumour increases sensitivity and "
    "shallower normal reduces the power to exclude shared variants."
)
CAPTION_C = (
    "Table 5.X. Clonal bottleneck diagnostic for the paired B4 call sets. "
    "Among variants at the Population tier, the number of alt-supporting "
    "reads observed in the matched wild-type sample. A variant with any "
    "read support in the wild-type may have pre-existed the clone; the "
    "Private tier retains only variants with zero support and is the set "
    "used for downstream analysis. Single-read support at ~30\u00D7 is "
    "within sequencing-error expectation and is not on its own evidence "
    "of pre-existence."
)


def load_cascade(cascade_dir: Path) -> dict:
    recs = {}
    for s in SAMPLES:
        p = cascade_dir / f"{s}_cascade.json"
        if not p.exists():
            sys.exit(f"Missing {p.relative_to(REPO)}; run "
                     "scripts/qc/filter_cascade.py --all first.")
        with p.open() as fh:
            recs[s] = json.load(fh)
    cs = {r.get("contig_set", "all") for r in recs.values()}
    if cs != {"primary"}:
        print(f"WARNING: cascade JSONs have contig_set={cs}; captions state "
              "primary contigs.")
    return recs


def load_depth(tsv: Path) -> dict:
    if not tsv.exists():
        sys.exit(f"Missing {tsv.relative_to(REPO)}; run "
                 "scripts/qc/bam_metrics.py --table first.")
    out = {}
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["sample"]] = float(row["mean_coverage"])
    return out


def final_tier(rec: dict) -> str:
    return "T5_PRIVATE" if rec["paired"] else "T4_POP"


def add_caption(doc, text):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.size = Pt(10)


def add_table(doc, header, rows, right_align_from=1, first_col_cm=3.0,
              total_cm=16.0):
    t = doc.add_table(rows=1 + len(rows), cols=len(header))
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.style = "Light Grid Accent 1"
    t.autofit = False
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    t._tbl.tblPr.append(layout)
    other = (total_cm - first_col_cm) / (len(header) - 1)
    for ci, col in enumerate(t.columns):
        col.width = Cm(first_col_cm if ci == 0 else other)
    for r in t.rows:
        for ci, c in enumerate(r.cells):
            c.width = Cm(first_col_cm if ci == 0 else other)
    for i, h in enumerate(header):
        t.rows[0].cells[i].text = h
        for p in t.rows[0].cells[i].paragraphs:
            for run in p.runs:
                run.bold = True
    for ri, row in enumerate(rows, start=1):
        for ci, val in enumerate(row):
            cell = t.rows[ri].cells[ci]
            cell.text = val
            if ci >= right_align_from:
                for p in cell.paragraphs:
                    p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    for r in t.rows:
        for c in r.cells:
            for p in c.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(9)
    doc.add_paragraph()


def fmt_int(v):
    return "\u2014" if v is None else f"{int(v):,}"


def fmt_tstv(v):
    return "\u2014" if v is None else f"{v:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cascade-dir", type=Path, default=CASCADE_DIR)
    ap.add_argument("--bam-tsv", type=Path, default=BAM_TSV)
    ap.add_argument("--out", type=Path, default=DOCX_OUT)
    a = ap.parse_args()

    recs = load_cascade(a.cascade_dir)
    depth = load_depth(a.bam_tsv)

    doc = Document()

    # ---- Table A: counts by tier ------------------------------------- #
    add_caption(doc, CAPTION_A)
    header = ["Sample"] + TIER_HEADERS + ["Ts/Tv (PASS)", "Ts/Tv (final)"]
    rows = []
    for s in SAMPLES:
        r = recs[s]
        ft = final_tier(r)
        cells = [LABELS[s]]
        for t in TIERS:
            n = r["tiers"][t]["n"]
            if t == "T5_PRIVATE" and not r["paired"]:
                cells.append("\u2014")
            else:
                cells.append(fmt_int(n))
        cells.append(fmt_tstv(r["tiers"]["T1_PASS"]["tstv"]))
        cells.append(fmt_tstv(r["tiers"][ft]["tstv"]))
        rows.append(cells)
    add_table(doc, header, rows)

    # ---- Table B: treatment contrast at final tier + depth ------------- #
    add_caption(doc, CAPTION_B)
    header = ["Genotype", "Tier", "Untreated", "Treated", "Change (%)",
              "Depth untreated (\u00D7)", "Depth treated (\u00D7)",
              "Tumour/normal depth ratio"]
    rows = []
    for geno, (un, tr, normals) in {
        "WT": ("WTUN", "WTAPH", None),
        "B4": ("B4UN", "B4APH", ("WTUN", "WTAPH")),
    }.items():
        ft = final_tier(recs[un])
        n_un = recs[un]["tiers"][ft]["n"]
        n_tr = recs[tr]["tiers"][ft]["n"]
        change = 100.0 * (n_tr - n_un) / n_un if n_un else float("nan")
        if normals:
            ratio = (f"{depth[un] / depth[normals[0]]:.2f} / "
                     f"{depth[tr] / depth[normals[1]]:.2f}")
        else:
            ratio = "\u2014 (tumour-only)"
        rows.append([
            geno, TIER_LABELS[TIERS.index(ft)],
            fmt_int(n_un), fmt_int(n_tr), f"{change:+.1f}",
            f"{depth[un]:.2f}", f"{depth[tr]:.2f}", ratio,
        ])
    add_table(doc, header, rows, right_align_from=2, first_col_cm=2.6)

    # ---- Table C: bottleneck --------------------------------------------- #
    add_caption(doc, CAPTION_C)
    header = ["Sample", "Population-tier variants", "0 reads", "1 read",
              "2 reads", "\u2265 3 reads", "Any support (%)"]
    rows = []
    for s in ("B4UN", "B4APH"):
        ns = recs[s].get("normal_alt_support_at_T4")
        if not ns:
            continue
        tot = sum(ns.values())
        sup = tot - ns["0"]
        rows.append([
            LABELS[s], fmt_int(tot), fmt_int(ns["0"]), fmt_int(ns["1"]),
            fmt_int(ns["2"]), fmt_int(ns["3+"]),
            f"{100.0 * sup / tot:.2f}" if tot else "\u2014",
        ])
    add_table(doc, header, rows)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(a.out)
    print(f"Wrote {a.out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
