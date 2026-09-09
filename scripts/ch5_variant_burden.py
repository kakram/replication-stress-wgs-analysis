"""
Chapter 5 §5.3 — Genome-wide somatic variant burden in MCF-7.

Thin wrapper around scripts/lib/somatic_burden.py. Defines the four MCF-7
sample paths and condition labels, invokes the shared analytical core,
applies MCF-7-specific output paths and table captions, and prints the
WT-vs-B4 ratio diagnostic the §5.3 prose is built on.

Outputs:
  outputs/ch5_variant_burden.csv
  outputs/ch5_variant_burden_thesis_table.tex
  outputs/ch5_variant_burden_summary.md

Usage:
  python3 scripts/ch5_variant_burden.py
"""

from __future__ import annotations

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

CSV_OUT = OUT_DIR / "ch5_variant_burden.csv"
TEX_OUT = OUT_DIR / "ch5_variant_burden_thesis_table.tex"
MD_OUT = OUT_DIR / "ch5_variant_burden_summary.md"

CAPTION = (
    "Table 5.3.1. Genome-wide somatic variant burden across four MCF-7 "
    "conditions. Total variants, PASS-filtered counts, filter-class "
    "breakdown, coding consequence distribution, SIFT predictions, "
    "gnomAD novelty, and tumour mutational burden (TMB) are shown for "
    "each sample. Filter labels follow GATK/Mutect2 conventions; rows "
    "carrying multiple filter terms increment each matching filter "
    "column once, so column sums may exceed (total $-$ PASS)."
)
LABEL = "tab:ch5-variant-burden"


def diagnostic_paragraph(results_by_name: dict[str, SampleResult]) -> str:
    """One-paragraph interpretation of the WT-vs-B4 asymmetry pattern,
    derived purely from the numbers we just computed."""
    wt_u, wt_a = results_by_name["WTUN"], results_by_name["WTAPH"]
    b4_u, b4_a = results_by_name["B4UN"], results_by_name["B4APH"]

    def ratio(a, b):
        return a / b if b else float("nan")

    r_un_pass = ratio(wt_u.pass_variants, b4_u.pass_variants)
    r_aph_pass = ratio(wt_a.pass_variants, b4_a.pass_variants)
    r_un_tmb = ratio(wt_u.tmb_mut_per_mb, b4_u.tmb_mut_per_mb)
    r_aph_tmb = ratio(wt_a.tmb_mut_per_mb, b4_a.tmb_mut_per_mb)

    # Compare the untreated and aphidicolin-treated ratios. If they are
    # within ±20 % of each other the asymmetry is treatment-independent;
    # otherwise APH modulates the WT-vs-B4 gap.
    def rel(a, b):
        if a is None or b is None or a != a or b != b:  # NaN check
            return float("inf")
        return abs(a - b) / ((a + b) / 2) if (a + b) else float("inf")

    consistent = rel(r_un_pass, r_aph_pass) <= 0.20

    direction_u = "WT > B4" if r_un_pass > 1 else "B4 > WT"
    direction_a = "WT > B4" if r_aph_pass > 1 else "B4 > WT"

    if consistent:
        verdict = (
            f"The WT-vs-B4 asymmetry is **consistent across treatment**: "
            f"the PASS-variant ratio is {r_un_pass:.2f} untreated vs "
            f"{r_aph_pass:.2f} aphidicolin-treated (within ±20 %), and "
            f"the TMB ratio shows the same pattern ({r_un_tmb:.2f} vs "
            f"{r_aph_tmb:.2f}). This is the signature of a **baseline-"
            f"burden / clonal-heterogeneity effect** rather than a "
            f"genotype × treatment interaction: the WT and B4 lines "
            f"differ in their somatic-variant load irrespective of "
            f"aphidicolin exposure. Plausible causes include the WT "
            f"line being a polyclonal MCF-7 stock against which the B4 "
            f"single-cell-derived ZFP36L1$^{{-/-}}$ clone is a "
            f"bottlenecked subline (so 'WT > B4' may simply reflect "
            f"pre-existing clonal diversity), or systematic differences "
            f"in coverage / variant-calling sensitivity between the two "
            f"line states. The §5.3 prose should foreground this "
            f"clonality caveat before claiming any aphidicolin or "
            f"ZFP36L1-loss effect from these counts."
        )
    else:
        verdict = (
            f"The WT-vs-B4 ratio **changes with aphidicolin treatment**: "
            f"PASS-variant ratio shifts from {r_un_pass:.2f} ({direction_u}) "
            f"untreated to {r_aph_pass:.2f} ({direction_a}) treated; TMB "
            f"ratio shifts from {r_un_tmb:.2f} to {r_aph_tmb:.2f}. This "
            f"pattern is consistent with a **genotype $\\times$ treatment "
            f"interaction**: the magnitude (and possibly direction) of "
            f"the WT-vs-B4 burden difference depends on whether the cells "
            f"have been exposed to replication stress, which is the "
            f"biologically informative outcome for the §5.3 narrative. "
            f"Note however that single replicates make any interaction "
            f"claim provisional — the §5.3 prose should flag this and "
            f"recommend technical replication before quantitative claims."
        )
    return verdict


def main() -> int:
    if not DATA_DIR.exists():
        raise SystemExit(f"FATAL: {DATA_DIR} not found. Copy MCF-7 sample "
                         "files from the external drive first.")

    print("=" * 72)
    print("Chapter 5 §5.3 — MCF-7 somatic variant burden")
    print("=" * 72)

    results: list[SampleResult] = []
    issues: list[str] = []
    for sample, _label in MCF7_SAMPLES:
        sample_dir = DATA_DIR / sample
        try:
            r = analyse_sample(sample, sample_dir, verbose=True)
            results.append(r)
        except FileNotFoundError as e:
            issues.append(str(e))
            print(f"[{sample}] SKIPPED: {e}")
            continue

    if len(results) != len(MCF7_SAMPLES):
        print()
        print("WARNING: some samples skipped; proceeding with what we have.")

    by_name = {r.sample: r for r in results}

    # CSV
    write_csv(results, CSV_OUT)
    print(f"\nCSV   -> {CSV_OUT.relative_to(REPO)}")

    # LaTeX
    write_latex_table(results, TEX_OUT, caption=CAPTION, label=LABEL)
    print(f"LaTeX -> {TEX_OUT.relative_to(REPO)}")

    # WT-vs-B4 ratio block + diagnostic
    ratio_lines: list[str] = []
    if {"WTUN", "B4UN"} <= by_name.keys():
        ratio_lines.append(
            pair_ratio_block(by_name["WTUN"], by_name["B4UN"],
                             "WT-U", "B4-U")
        )
    if {"WTAPH", "B4APH"} <= by_name.keys():
        ratio_lines.append(
            pair_ratio_block(by_name["WTAPH"], by_name["B4APH"],
                             "WT-A", "B4-A")
        )

    print()
    print("WT-vs-B4 ratio diagnostic")
    print("-" * 72)
    for block in ratio_lines:
        print(block)
        print()

    per_sample_summary = (
        f"{'Sample':<8}  {'PASS':>10}  {'Coding':>8}  {'Missense':>9}  "
        f"{'SIFT-del':>9}  {'TMB':>6}\n"
    )
    for r in results:
        per_sample_summary += (
            f"{r.sample:<8}  {r.pass_variants:>10,}  "
            f"{r.coding_total:>8,}  {r.missense:>9,}  "
            f"{r.sift_deleterious:>9,}  {r.tmb_mut_per_mb:>6.2f}\n"
        )
    print("Per-sample summary")
    print("-" * 72)
    print(per_sample_summary)

    diagnostic = diagnostic_paragraph(by_name) if len(by_name) == 4 else (
        "Insufficient samples to compute the WT-vs-B4 × treatment "
        "diagnostic (need all four samples)."
    )
    print("Diagnostic")
    print("-" * 72)
    print(diagnostic)

    # Markdown summary
    md_lines: list[str] = []
    md_lines.append("# §5.3 MCF-7 Somatic Variant Burden — Summary\n")
    md_lines.append(
        "Source data: Sentieon TNhaplotyper2 + TNfilter somatic calls for "
        "MCF-7 WT (parental) and B4 (CRISPR/Cas9 ZFP36L1$^{-/-}$ clone, "
        "from Teotia 2024), each $\\pm$ aphidicolin. Files used per "
        "sample: `<SAMPLE>_somatic.vcf.gz` (FILTER classes), "
        "`<SAMPLE>_somatic_vep_anno.maf.gz` (coding / SIFT / gnomAD), "
        "`<SAMPLE>_tmb.tsv` (TMB).\n"
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
    MD_OUT.write_text("\n".join(md_lines), encoding="utf-8")
    print(f"\nMD    -> {MD_OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
