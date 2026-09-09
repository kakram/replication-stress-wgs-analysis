"""
Chapter 5 §5.3 — Genome-wide somatic variant burden in MCF-7.

Thin wrapper around scripts/lib/somatic_burden.py. Defines the four MCF-7
sample paths and condition labels, invokes the shared analytical core,
applies MCF-7-specific output paths and table captions, and prints the
WT-vs-B4 ratio diagnostic the §5.3 prose is built on.

Contig policy (2026-09-09): counts are restricted to the primary assembly
(chr1-22, chrX, chrY, chrM) by default. Pass --all-contigs to reproduce
the original all-contig numbers; outputs then carry an `_allcontigs`
suffix so the two never overwrite each other.

Outputs (default, primary contigs):
  outputs/ch5_variant_burden.csv
  outputs/ch5_variant_burden_thesis_table.tex
  outputs/ch5_variant_burden_summary.md

Usage:
  python3 scripts/ch5_variant_burden.py [--all-contigs]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from lib.somatic_burden import (  # noqa: E402
    SampleResult,
    analyse_sample,
    markdown_table,
    pair_ratio_block,
    write_csv,
    write_latex_table,
)


# MCF-7 sample definitions. Order is meaningful: ratios are computed within
# each (WT, B4) pair under matched aphidicolin status.
MCF7_SAMPLES = [
    ("WTUN",  "WT untreated"),
    ("WTAPH", "WT + aphidicolin"),
    ("B4UN",  "B4 ZFP36L1$^{-/-}$ untreated"),
    ("B4APH", "B4 ZFP36L1$^{-/-}$ + aphidicolin"),
]

DATA_DIR = REPO / "data" / "mcf7"
OUT_DIR = REPO / "outputs"

CAPTION = (
    "Table 5.X. Genome-wide somatic small-variant burden across four MCF-7 "
    "conditions. Wild-type samples were called in tumour-only mode and B4 "
    "samples in paired tumour--normal mode against the treatment-matched "
    "wild-type, so counts are comparable within, but not between, "
    "genotypes. Counts are restricted to the primary assembly (chr1--22, "
    "X, Y, M); records on unplaced, alternate, decoy and HLA contigs are "
    "excluded and their number is shown. Total variants, PASS-filtered "
    "counts, filter-class breakdown, coding consequence distribution, SIFT "
    "predictions and gnomAD novelty are shown for each sample. Filter "
    "labels follow GATK/Mutect2 conventions; rows carrying multiple filter "
    "terms increment each matching filter column once, so column sums may "
    "exceed (total $-$ PASS)."
)
LABEL = "tab:ch5-variant-burden"


def diagnostic_paragraph(results_by_name: dict[str, SampleResult]) -> str:
    """One-paragraph interpretation of the WT-vs-B4 pattern, derived from
    the numbers just computed. Rewritten 2026-09-09 to the calling-design
    explanation established in the 2026-09-02 verification report; the
    earlier clonal-heterogeneity wording is superseded."""
    wt_u, wt_a = results_by_name["WTUN"], results_by_name["WTAPH"]
    b4_u, b4_a = results_by_name["B4UN"], results_by_name["B4APH"]

    def ratio(a, b):
        return a / b if b else float("nan")

    r_un = ratio(wt_u.pass_variants, b4_u.pass_variants)
    r_aph = ratio(wt_a.pass_variants, b4_a.pass_variants)
    d_wt = 100.0 * (wt_a.pass_variants - wt_u.pass_variants) / wt_u.pass_variants
    d_b4 = 100.0 * (b4_a.pass_variants - b4_u.pass_variants) / b4_u.pass_variants

    def rel(a, b):
        return abs(a - b) / ((a + b) / 2) if (a + b) else float("inf")

    invariant = rel(r_un, r_aph) <= 0.20

    lines = [
        f"WT-to-B4 PASS ratio: {r_un:.2f} untreated, {r_aph:.2f} "
        f"aphidicolin-treated"
        + (" (within 20% of each other: treatment-invariant)." if invariant
           else " (differ by more than 20%: check before interpreting)."),
        "",
        "The gap is a property of the calling design, not of the biology. "
        "B4 sets were called in paired mode against the treatment-matched "
        "wild-type and contain only variants private to the clone; wild-type "
        "sets were called tumour-only with no matched normal and no panel of "
        "normals, and retain the cell line's full complement of variation "
        "not modelled as germline by gnomAD. The polyclonal architecture of "
        "the wild-type stock and the monoclonal origin of B4 are additional "
        "contributors that these data cannot separate. WT-vs-B4 burden must "
        "not be interpreted as an effect of ZFP36L1 loss.",
        "",
        f"Within-genotype treatment contrast (the only valid burden "
        f"comparison here): WT {wt_u.pass_variants:,} -> "
        f"{wt_a.pass_variants:,} ({d_wt:+.1f}%); B4 {b4_u.pass_variants:,} "
        f"-> {b4_a.pass_variants:,} ({d_b4:+.1f}%). B4UN was sequenced at "
        f"44.78x and B4APH at 34.33x; see the depth-matched filter cascade "
        f"(results/cascade_primary) before reading any B4 change as a "
        f"treatment effect.",
    ]
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-contigs", action="store_true",
                    help="count every contig (original behaviour); outputs "
                         "get an _allcontigs suffix")
    args = ap.parse_args()
    primary_only = not args.all_contigs
    suffix = "" if primary_only else "_allcontigs"

    csv_out = OUT_DIR / f"ch5_variant_burden{suffix}.csv"
    tex_out = OUT_DIR / f"ch5_variant_burden_thesis_table{suffix}.tex"
    md_out = OUT_DIR / f"ch5_variant_burden_summary{suffix}.md"

    if not DATA_DIR.exists():
        raise SystemExit(f"FATAL: {DATA_DIR} not found. Copy MCF-7 sample "
                         "files from the external drive first.")

    print("=" * 72)
    print("Chapter 5 §5.3 — MCF-7 somatic variant burden "
          f"[{'primary contigs' if primary_only else 'ALL contigs'}]")
    print("=" * 72)

    results: list[SampleResult] = []
    issues: list[str] = []
    for sample, _label in MCF7_SAMPLES:
        sample_dir = DATA_DIR / sample
        try:
            r = analyse_sample(sample, sample_dir, verbose=True,
                               primary_only=primary_only)
            results.append(r)
        except FileNotFoundError as e:
            issues.append(str(e))
            print(f"[{sample}] SKIPPED: {e}")
            continue

    if len(results) != len(MCF7_SAMPLES):
        print()
        print("WARNING: some samples skipped; proceeding with what we have.")

    by_name = {r.sample: r for r in results}

    write_csv(results, csv_out)
    print(f"\nCSV   -> {csv_out.relative_to(REPO)}")

    write_latex_table(results, tex_out, caption=CAPTION, label=LABEL)
    print(f"LaTeX -> {tex_out.relative_to(REPO)}")

    ratio_lines: list[str] = []
    if {"WTUN", "B4UN"} <= by_name.keys():
        ratio_lines.append(
            pair_ratio_block(by_name["WTUN"], by_name["B4UN"], "WT-U", "B4-U")
        )
    if {"WTAPH", "B4APH"} <= by_name.keys():
        ratio_lines.append(
            pair_ratio_block(by_name["WTAPH"], by_name["B4APH"], "WT-A", "B4-A")
        )

    print()
    print("WT-vs-B4 ratio diagnostic")
    print("-" * 72)
    for block in ratio_lines:
        print(block)
        print()

    per_sample_summary = (
        f"{'Sample':<8}  {'Total':>10}  {'PASS':>10}  {'Excl.':>8}  "
        f"{'Coding':>8}  {'Missense':>9}  {'SIFT-del':>9}\n"
    )
    for r in results:
        per_sample_summary += (
            f"{r.sample:<8}  {r.total_variants:>10,}  {r.pass_variants:>10,}  "
            f"{r.vcf_records_nonprimary:>8,}  {r.coding_total:>8,}  "
            f"{r.missense:>9,}  {r.sift_deleterious:>9,}\n"
        )
    print("Per-sample summary (Excl. = records on non-primary contigs)")
    print("-" * 72)
    print(per_sample_summary)

    diagnostic = diagnostic_paragraph(by_name) if len(by_name) == 4 else (
        "Insufficient samples to compute the WT-vs-B4 diagnostic "
        "(need all four samples)."
    )
    print("Diagnostic")
    print("-" * 72)
    print(diagnostic)

    md_lines: list[str] = []
    md_lines.append("# §5.3 MCF-7 Somatic Variant Burden — Summary\n")
    md_lines.append(
        f"Contig set: **{'primary (chr1-22, X, Y, M)' if primary_only else 'all'}**. "
        "Source data: Sentieon TNhaplotyper2 + TNfilter somatic calls for "
        "MCF-7 WT (parental; tumour-only calling) and B4 (CRISPR/Cas9 "
        "ZFP36L1$^{-/-}$ clone from Teotia 2024; paired calling against the "
        "treatment-matched WT), each $\\pm$ aphidicolin. Files used per "
        "sample: `<SAMPLE>_somatic.vcf.gz` (FILTER classes), "
        "`<SAMPLE>_somatic_vep_anno.maf.gz` (coding / SIFT / gnomAD). "
        "Provider TMB is carried in the CSV but not reported.\n"
    )
    md_lines.append("## Numbers\n")
    md_lines.append(markdown_table(results) + "\n")
    md_lines.append("## WT-vs-B4 ratios\n")
    md_lines.append("```\n" + "\n\n".join(ratio_lines) + "\n```\n")
    md_lines.append("## Per-sample summary\n")
    md_lines.append("```\n" + per_sample_summary + "```\n")
    md_lines.append("## Diagnostic\n")
    md_lines.append(diagnostic + "\n")
    if issues:
        md_lines.append("## Data issues\n")
        for it in issues:
            md_lines.append(f"- {it}\n")
    md_out.write_text("\n".join(md_lines), encoding="utf-8")
    print(f"\nMD    -> {md_out.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
