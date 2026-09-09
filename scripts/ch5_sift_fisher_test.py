#!/usr/bin/env python3
"""
SIFT del:tol Fisher exact testing — §5.10 (referenced from §5.3.1).

§5.3.1 paragraph 3 identified a divergence in the ratio of SIFT-deleterious
to SIFT-tolerated missense calls between the MCF-7 wild-type (WT) parental
stock (~1 : 1.9 in both untreated and aphidicolin-treated samples) and
the B4 ZFP36L1-/- clone (~1 : 1 in both). With only ~65 SIFT-callable
missense calls per B4 condition, the §5.3.1 prose deferred formal
testing to the present (§5.10) analysis.

This script implements a hierarchical Fisher's-exact testing strategy:

  Test 1: WTUN vs WTAPH       — within-WT treatment effect
  Test 2: B4UN vs B4APH       — within-B4 treatment effect
  Test 3: WT pooled vs B4 pooled — between-genotype effect
  Test 4: 2x4 chi-square      — global test across all four samples

Logic: Tests 1 and 2 establish whether treatment is a confound. If
neither rejects (as expected from §5.3.1 paragraph 3), pooling across
treatment status for Test 3 is justified. Test 4 is a global sanity
check that the pattern from Tests 1-3 holds in a single omnibus
contrast. Bonferroni correction is applied across the three pairwise
tests (alpha_corr = 0.05 / 3 = 0.0167); Test 4 is reported uncorrected
as a confirmatory check.

Effect sizes are reported as odds ratios with 95 % confidence intervals,
computed using scipy.stats.contingency.odds_ratio where available
(scipy >= 1.10), with a Wald log-odds-ratio fallback for older scipy
installations.

Inputs:
  outputs/ch5_variant_burden.csv

Outputs:
  outputs/ch5_sift_fisher_results.csv   — machine-readable
  outputs/ch5_sift_fisher_summary.md    — markdown for §5.10 drafting

Usage:
    python3 scripts/ch5_sift_fisher_test.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact


REPO = Path(__file__).resolve().parent.parent
CSV_IN = REPO / "outputs" / "ch5_variant_burden.csv"
CSV_OUT = REPO / "outputs" / "ch5_sift_fisher_results.csv"
MD_OUT = REPO / "outputs" / "ch5_sift_fisher_summary.md"

SAMPLES = ["WTUN", "WTAPH", "B4UN", "B4APH"]


def or_with_ci(table):
    """Return (odds_ratio, ci_low, ci_high) for a 2x2 contingency table.

    Prefers scipy.stats.contingency.odds_ratio (exact Fisher CI, scipy >=
    1.10). Falls back to the Wald log-odds-ratio CI for older scipy.
    """
    try:
        from scipy.stats.contingency import odds_ratio as scipy_or
        result = scipy_or(table)
        ci = result.confidence_interval(confidence_level=0.95)
        return float(result.statistic), float(ci.low), float(ci.high)
    except (ImportError, AttributeError):
        a, b = table[0]
        c, d = table[1]
        if min(a, b, c, d) == 0:
            return float("nan"), float("nan"), float("nan")
        or_val = (a * d) / (b * c)
        se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        log_or = np.log(or_val)
        return or_val, float(np.exp(log_or - 1.96 * se)), float(np.exp(log_or + 1.96 * se))


def fisher_block(label: str, table):
    """Run Fisher's exact test on a 2x2 table and assemble a results row."""
    a, b = table[0]
    c, d = table[1]
    _, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    or_val, ci_low, ci_high = or_with_ci([[a, b], [c, d]])
    return {
        "test": label,
        "g1_del": a, "g1_tol": b,
        "g2_del": c, "g2_tol": d,
        "g1_pct_del": 100 * a / (a + b),
        "g2_pct_del": 100 * c / (c + d),
        "odds_ratio": or_val,
        "ci95_low": ci_low,
        "ci95_high": ci_high,
        "p_value": p_value,
    }


def main():
    if not CSV_IN.exists():
        sys.exit(f"Input CSV not found: {CSV_IN}")

    df = pd.read_csv(CSV_IN).set_index("sample")

    counts = {
        s: {"del": int(df.loc[s, "sift_deleterious"]),
            "tol": int(df.loc[s, "sift_tolerated"])}
        for s in SAMPLES
    }
    wt_pool = {
        "del": counts["WTUN"]["del"] + counts["WTAPH"]["del"],
        "tol": counts["WTUN"]["tol"] + counts["WTAPH"]["tol"],
    }
    b4_pool = {
        "del": counts["B4UN"]["del"] + counts["B4APH"]["del"],
        "tol": counts["B4UN"]["tol"] + counts["B4APH"]["tol"],
    }

    # ── Console output ────────────────────────────────────────────────────
    print("=== SIFT del:tol counts per sample ===\n")
    for s in SAMPLES:
        n = counts[s]["del"] + counts[s]["tol"]
        pct = 100 * counts[s]["del"] / n
        print(f"  {s:6s}  del={counts[s]['del']:>4}  tol={counts[s]['tol']:>4}  "
              f"total={n:>4}  pct_del={pct:5.1f}%")
    wt_total = wt_pool["del"] + wt_pool["tol"]
    b4_total = b4_pool["del"] + b4_pool["tol"]
    print(f"\n  WT pool del={wt_pool['del']:>4}  tol={wt_pool['tol']:>4}  "
          f"total={wt_total:>4}  pct_del={100*wt_pool['del']/wt_total:5.1f}%")
    print(f"  B4 pool del={b4_pool['del']:>4}  tol={b4_pool['tol']:>4}  "
          f"total={b4_total:>4}  pct_del={100*b4_pool['del']/b4_total:5.1f}%")

    alpha_corr = 0.05 / 3
    print(f"\n=== Fisher's exact tests (Bonferroni \u03B1_corr = {alpha_corr:.4f}) ===\n")

    results = [
        fisher_block(
            "Test 1: WTUN vs WTAPH (within-WT treatment)",
            [[counts["WTUN"]["del"], counts["WTUN"]["tol"]],
             [counts["WTAPH"]["del"], counts["WTAPH"]["tol"]]],
        ),
        fisher_block(
            "Test 2: B4UN vs B4APH (within-B4 treatment)",
            [[counts["B4UN"]["del"], counts["B4UN"]["tol"]],
             [counts["B4APH"]["del"], counts["B4APH"]["tol"]]],
        ),
        fisher_block(
            "Test 3: WT pooled vs B4 pooled (between-genotype)",
            [[wt_pool["del"], wt_pool["tol"]],
             [b4_pool["del"], b4_pool["tol"]]],
        ),
    ]

    for r in results:
        sig = "SIGNIFICANT" if r["p_value"] < alpha_corr else "not significant"
        print(f"  {r['test']}")
        print(f"    Group 1: {r['g1_del']}/{r['g1_del']+r['g1_tol']} = "
              f"{r['g1_pct_del']:5.1f}% deleterious")
        print(f"    Group 2: {r['g2_del']}/{r['g2_del']+r['g2_tol']} = "
              f"{r['g2_pct_del']:5.1f}% deleterious")
        print(f"    OR = {r['odds_ratio']:.3f}  "
              f"(95% CI: {r['ci95_low']:.3f}, {r['ci95_high']:.3f})")
        print(f"    p  = {r['p_value']:.4g}  -> {sig}\n")

    # ── Test 4: 2x4 chi-square ────────────────────────────────────────────
    table_2x4 = np.array([
        [counts["WTUN"]["del"], counts["WTAPH"]["del"],
         counts["B4UN"]["del"], counts["B4APH"]["del"]],
        [counts["WTUN"]["tol"], counts["WTAPH"]["tol"],
         counts["B4UN"]["tol"], counts["B4APH"]["tol"]],
    ])
    chi2, p_chi2, dof, expected = chi2_contingency(table_2x4)
    print(f"  Test 4: 2x4 chi-square (sample x SIFT class)")
    print(f"    chi2 = {chi2:.3f},  dof = {dof},  p = {p_chi2:.4g}\n")

    # ── CSV output ───────────────────────────────────────────────────────
    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(CSV_OUT, index=False)
    print(f"Wrote {CSV_OUT.relative_to(REPO)}")

    # ── Markdown summary ─────────────────────────────────────────────────
    md = [
        "# \u00A75.10 SIFT del:tol Fisher exact tests \u2014 Summary",
        "",
        "Source: `outputs/ch5_variant_burden.csv` (sift_deleterious, sift_tolerated).",
        "",
        "## SIFT counts per sample",
        "",
        "| Sample | SIFT-del | SIFT-tol | Total | % del |",
        "|---|---:|---:|---:|---:|",
    ]
    for s in SAMPLES:
        n = counts[s]["del"] + counts[s]["tol"]
        pct = 100 * counts[s]["del"] / n
        md.append(f"| {s} | {counts[s]['del']} | {counts[s]['tol']} | {n} | {pct:.1f}% |")
    md.append(f"| **WT pooled** | **{wt_pool['del']}** | **{wt_pool['tol']}** | "
              f"**{wt_total}** | **{100*wt_pool['del']/wt_total:.1f}%** |")
    md.append(f"| **B4 pooled** | **{b4_pool['del']}** | **{b4_pool['tol']}** | "
              f"**{b4_total}** | **{100*b4_pool['del']/b4_total:.1f}%** |")
    md.append("")
    md.append(f"Bonferroni-corrected threshold for the three pairwise tests: "
              f"\u03B1\u209C\u2092\u1D63\u1D63 = 0.05 / 3 = **{alpha_corr:.4f}**.")
    md.append("")
    md.append("## Test results")
    md.append("")
    for r in results:
        sig = "**SIGNIFICANT**" if r["p_value"] < alpha_corr else "not significant"
        md.append(f"### {r['test']}")
        md.append("")
        md.append(f"- Group 1: {r['g1_del']} / {r['g1_del']+r['g1_tol']} "
                  f"({r['g1_pct_del']:.1f}%) deleterious")
        md.append(f"- Group 2: {r['g2_del']} / {r['g2_del']+r['g2_tol']} "
                  f"({r['g2_pct_del']:.1f}%) deleterious")
        md.append(f"- Odds ratio: **{r['odds_ratio']:.3f}**  "
                  f"(95% CI: {r['ci95_low']:.3f}, {r['ci95_high']:.3f})")
        md.append(f"- p-value: **{r['p_value']:.4g}** \u2014 {sig}")
        md.append("")
    md.append("### Test 4: 2x4 chi-square across all four samples")
    md.append("")
    md.append(f"- \u03C7\u00B2 = {chi2:.3f}, dof = {dof}, p = {p_chi2:.4g}")
    md.append("")

    MD_OUT.write_text("\n".join(md))
    print(f"Wrote {MD_OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
