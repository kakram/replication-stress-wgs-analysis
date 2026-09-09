"""
Somatic variant burden analysis library.

Shared analytical core for §5.3 (MCF-7, Sentieon TNhaplotyper2 + TNfilter)
and any future somatic re-run of Chapter 6 (U2OS). Kept free of
cell-line-specific assumptions so the wrappers stay thin.

Inputs per sample:
  - <SAMPLE>_somatic.vcf.gz      : source for FILTER-class counts.
                                    MAF FILTER column collapses to PASS /
                                    common_variant; the VCF preserves the
                                    full Mutect2 filter vocabulary
                                    (germline, panel_of_normals,
                                    weak_evidence, clustered_events, etc.).
  - <SAMPLE>_somatic_vep_anno.maf.gz : source for coding consequences,
                                    SIFT, gnomAD novelty, gene counts.
  - <SAMPLE>_tmb.tsv             : pre-computed TMB (mut/Mb).

Filter-counting policy: a multi-filter VCF row (e.g.
'clustered_events;germline;weak_evidence') increments **each** matching
named-filter bucket once, plus filter_other if any token isn't one of the
four named filters. Consequence: column sums may exceed (total - PASS)
when rows have multiple filters, which is the desired epidemiological
behaviour.
"""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from dataclasses import dataclass, field, fields as dc_fields
from pathlib import Path
from typing import Iterable, Sequence


# MAF Variant_Classification values we treat as "coding".
CODING_CLASSES = {
    "Missense_Mutation",
    "Silent",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Splice_Region",
    "Nonstop_Mutation",
    "Translation_Start_Site",
}

CLASS_TO_FIELD = {
    "Missense_Mutation":   "missense",
    "Silent":              "synonymous",
    "Nonsense_Mutation":   "nonsense",
    "Frame_Shift_Del":     "frameshift_del",
    "Frame_Shift_Ins":     "frameshift_ins",
    "In_Frame_Del":        "in_frame_del",
    "In_Frame_Ins":        "in_frame_ins",
    "Splice_Site":         "splice_site",
    "Splice_Region":       "splice_region",
}

NAMED_FILTER_FIELDS = {
    "germline":         "filter_germline",
    "panel_of_normals": "filter_panel_of_normals",
    "weak_evidence":    "filter_weak_evidence",
    "clustered_events": "filter_clustered",
}

# Threshold for "novel against gnomAD". gnomAD_AF strictly above this is
# treated as "known"; missing or <=threshold is "novel".
GNOMAD_NOVEL_THRESHOLD = 0.001


@dataclass
class SampleResult:
    sample: str
    # Whole-VCF counts
    total_variants: int = 0
    pass_variants: int = 0
    filter_germline: int = 0
    filter_panel_of_normals: int = 0
    filter_weak_evidence: int = 0
    filter_clustered: int = 0
    filter_other: int = 0
    # MAF coding consequences
    coding_total: int = 0
    missense: int = 0
    synonymous: int = 0
    nonsense: int = 0
    frameshift_del: int = 0
    frameshift_ins: int = 0
    in_frame_del: int = 0
    in_frame_ins: int = 0
    splice_site: int = 0
    splice_region: int = 0
    # Functional + novelty
    sift_deleterious: int = 0
    sift_tolerated: int = 0
    sift_unknown: int = 0
    gnomad_known: int = 0
    gnomad_novel: int = 0
    # External
    tmb_mut_per_mb: float = 0.0
    # Gene reach
    genes_with_any_variant: int = 0
    genes_with_coding_variant: int = 0


CSV_COLUMNS: Sequence[str] = tuple(f.name for f in dc_fields(SampleResult))


# ---------------------------- parsers ------------------------------------ #

def parse_vcf_filters(vcf_path: Path) -> dict:
    """Stream a bgzipped/gzipped VCF and tally filter classes."""
    counts = {
        "total_variants": 0,
        "pass_variants": 0,
        "filter_germline": 0,
        "filter_panel_of_normals": 0,
        "filter_weak_evidence": 0,
        "filter_clustered": 0,
        "filter_other": 0,
    }
    with gzip.open(str(vcf_path), "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            # FILTER is column 7 (1-indexed).
            parts = line.split("\t", 8)
            if len(parts) < 7:
                continue
            counts["total_variants"] += 1
            filt = parts[6]
            if filt == "PASS":
                counts["pass_variants"] += 1
                continue
            saw_named = False
            saw_unnamed = False
            for tok in filt.split(";"):
                field_name = NAMED_FILTER_FIELDS.get(tok)
                if field_name is not None:
                    counts[field_name] += 1
                    saw_named = True
                else:
                    saw_unnamed = True
            if saw_unnamed and not saw_named:
                counts["filter_other"] += 1
            elif saw_unnamed and saw_named:
                # Multi-filter row with at least one unnamed token: also
                # increment 'other' so the user sees that this row has
                # filters beyond the named four. Sum can exceed
                # (total - PASS) — documented behaviour.
                counts["filter_other"] += 1
    return counts


def parse_maf(maf_path: Path) -> dict:
    """Stream a MAF and tally consequence classes, SIFT, gnomAD, gene sets."""
    classes: Counter = Counter()
    sift_d = sift_t = sift_u = 0
    gnomad_known = gnomad_novel = 0
    genes_any: set[str] = set()
    genes_coding: set[str] = set()

    with gzip.open(str(maf_path), "rt") as fh:
        idx: dict[str, int] = {}
        header_seen = False
        for raw in fh:
            if not header_seen:
                if raw.startswith("#"):
                    continue
                header = raw.rstrip("\n").split("\t")
                idx = {name: i for i, name in enumerate(header)}
                header_seen = True
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < len(idx):
                fields += [""] * (len(idx) - len(fields))

            sym = fields[idx["Hugo_Symbol"]] if "Hugo_Symbol" in idx else ""
            cls = fields[idx["Variant_Classification"]] \
                if "Variant_Classification" in idx else ""
            sift = fields[idx["SIFT"]] if "SIFT" in idx else ""
            gnomad = fields[idx["gnomAD_AF"]] if "gnomAD_AF" in idx else ""

            classes[cls] += 1

            if sym and sym not in (".", ""):
                genes_any.add(sym)
                if cls in CODING_CLASSES:
                    genes_coding.add(sym)

            if sift.startswith("deleterious"):
                sift_d += 1
            elif sift.startswith("tolerated"):
                sift_t += 1
            else:
                sift_u += 1

            try:
                af = float(gnomad) if gnomad else 0.0
            except ValueError:
                af = 0.0
            if af > GNOMAD_NOVEL_THRESHOLD:
                gnomad_known += 1
            else:
                gnomad_novel += 1

    out = {
        "coding_total": sum(classes[c] for c in CODING_CLASSES if c in classes),
        "sift_deleterious": sift_d,
        "sift_tolerated": sift_t,
        "sift_unknown": sift_u,
        "gnomad_known": gnomad_known,
        "gnomad_novel": gnomad_novel,
        "genes_with_any_variant": len(genes_any),
        "genes_with_coding_variant": len(genes_coding),
    }
    for cls_name, fld in CLASS_TO_FIELD.items():
        out[fld] = classes.get(cls_name, 0)
    return out


def parse_tmb(tmb_path: Path) -> float:
    """Read the Sentieon TMB tsv (4 cols: sample, target_size, count, TMB).
    Returns TMB in mutations per Mb."""
    with tmb_path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return 0.0
    return float(rows[0]["TMB"])


# ---------------------------- per-sample driver -------------------------- #

def analyse_sample(sample_name: str, sample_dir: Path,
                   verbose: bool = True) -> SampleResult:
    """Run all three parsers for one sample and return a SampleResult.

    `sample_dir` is the per-sample directory containing the three expected
    files named <SAMPLE>_somatic.vcf.gz, <SAMPLE>_somatic_vep_anno.maf.gz
    and <SAMPLE>_tmb.tsv.
    """
    sample_dir = Path(sample_dir)
    vcf = sample_dir / f"{sample_name}_somatic.vcf.gz"
    maf = sample_dir / f"{sample_name}_somatic_vep_anno.maf.gz"
    tmb = sample_dir / f"{sample_name}_tmb.tsv"
    for f in (vcf, maf, tmb):
        if not f.exists():
            raise FileNotFoundError(
                f"Missing required file for sample {sample_name}: {f}"
            )

    result = SampleResult(sample=sample_name)

    if verbose:
        print(f"[{sample_name}] parsing VCF for filter classes "
              f"({vcf.stat().st_size / 1024**2:.1f} MB) …")
    for k, v in parse_vcf_filters(vcf).items():
        setattr(result, k, v)

    if verbose:
        print(f"[{sample_name}] parsing MAF for coding / SIFT / gnomAD "
              f"({maf.stat().st_size / 1024**2:.1f} MB) …")
    for k, v in parse_maf(maf).items():
        setattr(result, k, v)

    result.tmb_mut_per_mb = parse_tmb(tmb)

    if verbose:
        print(f"[{sample_name}] done: PASS={result.pass_variants:,}  "
              f"coding={result.coding_total:,}  "
              f"TMB={result.tmb_mut_per_mb}")
    return result


# ---------------------------- output writers ----------------------------- #

def write_csv(results: Iterable[SampleResult], csv_path: Path) -> None:
    csv_path = Path(csv_path)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for r in results:
            w.writerow([getattr(r, c) for c in CSV_COLUMNS])


# Pretty labels for the human-readable tables (LaTeX / markdown).
PRETTY_ROW_ORDER: Sequence[tuple[str, str]] = (
    ("total_variants",            r"Total variants (VCF)"),
    ("pass_variants",             r"PASS variants"),
    ("filter_germline",           r"\,\,FILTER: germline"),
    ("filter_panel_of_normals",   r"\,\,FILTER: panel\_of\_normals"),
    ("filter_weak_evidence",      r"\,\,FILTER: weak\_evidence"),
    ("filter_clustered",          r"\,\,FILTER: clustered\_events"),
    ("filter_other",              r"\,\,FILTER: other"),
    ("coding_total",              r"Coding variants (total)"),
    ("missense",                  r"\,\,Missense"),
    ("synonymous",                r"\,\,Synonymous"),
    ("nonsense",                  r"\,\,Nonsense"),
    ("frameshift_del",            r"\,\,Frameshift deletion"),
    ("frameshift_ins",            r"\,\,Frameshift insertion"),
    ("in_frame_del",              r"\,\,In-frame deletion"),
    ("in_frame_ins",              r"\,\,In-frame insertion"),
    ("splice_site",               r"\,\,Splice site"),
    ("splice_region",             r"\,\,Splice region"),
    ("sift_deleterious",          r"SIFT deleterious"),
    ("sift_tolerated",            r"SIFT tolerated"),
    ("sift_unknown",              r"SIFT unknown"),
    ("gnomad_known",              r"gnomAD known (AF > 0.001)"),
    ("gnomad_novel",              r"gnomAD novel (AF $\le$ 0.001)"),
    ("tmb_mut_per_mb",            r"TMB (mut/Mb)"),
    ("genes_with_any_variant",    r"Genes with $\ge$1 variant"),
    ("genes_with_coding_variant", r"Genes with coding variant"),
)


def write_latex_table(results: Sequence[SampleResult], tex_path: Path,
                      caption: str, label: str) -> None:
    """Emit a thesis-ready transposed table: metrics as rows, samples as
    columns. Used in §5.3 and re-usable for the future U2OS somatic run."""
    tex_path = Path(tex_path)
    tex_path.parent.mkdir(parents=True, exist_ok=True)
    samples = [r.sample for r in results]
    by_sample = {r.sample: r for r in results}

    col_spec = "l" + "r" * len(samples)
    lines: list[str] = []
    lines.append(r"\begin{table}[ht]")
    lines.append(r"\centering")
    lines.append(r"\small")
    lines.append(rf"\caption{{{caption}}}")
    lines.append(rf"\label{{{label}}}")
    lines.append(rf"\begin{{tabular}}{{{col_spec}}}")
    lines.append(r"\hline")
    header_cells = ["Metric"] + [s.replace("_", r"\_") for s in samples]
    lines.append(" & ".join(header_cells) + r" \\")
    lines.append(r"\hline")
    for field_name, pretty in PRETTY_ROW_ORDER:
        row_cells = [pretty]
        for s in samples:
            v = getattr(by_sample[s], field_name)
            if isinstance(v, float):
                row_cells.append(f"{v:.2f}")
            else:
                row_cells.append(f"{v:,}")
        lines.append(" & ".join(row_cells) + r" \\")
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    tex_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def markdown_table(results: Sequence[SampleResult]) -> str:
    """Return a markdown table mirroring the LaTeX one (no LaTeX escapes)."""
    samples = [r.sample for r in results]
    by_sample = {r.sample: r for r in results}
    lines = ["| Metric | " + " | ".join(samples) + " |"]
    lines.append("|" + "---|" * (len(samples) + 1))

    # Markdown-friendly labels: strip LaTeX dashes and replace indentation.
    def md_label(field_name: str, pretty_latex: str) -> str:
        s = pretty_latex.replace(r"\,\,", "  ")
        s = s.replace(r"\_", "_")
        s = s.replace(r"\le", "≤")
        s = s.replace(r"$\ge$", "≥").replace(r"\ge", "≥")
        s = s.replace("$", "")
        return s

    for field_name, pretty in PRETTY_ROW_ORDER:
        label = md_label(field_name, pretty)
        row = [label]
        for s in samples:
            v = getattr(by_sample[s], field_name)
            if isinstance(v, float):
                row.append(f"{v:.2f}")
            else:
                row.append(f"{v:,}")
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


# ---------------------------- ratio helpers ------------------------------ #

def _safe_ratio(a: float, b: float) -> str:
    if b == 0:
        return "n/a"
    return f"{a / b:.2f}"


def pair_ratio_block(a: SampleResult, b: SampleResult,
                     a_label: str, b_label: str) -> str:
    """Return a 4-line block of headline ratios between two samples."""
    return (
        f"{a_label} vs {b_label}:\n"
        f"  PASS variants ratio   : {a.pass_variants:>8,} / "
        f"{b.pass_variants:>8,}  = {_safe_ratio(a.pass_variants, b.pass_variants)}\n"
        f"  Coding variants ratio : {a.coding_total:>8,} / "
        f"{b.coding_total:>8,}  = {_safe_ratio(a.coding_total, b.coding_total)}\n"
        f"  TMB ratio             : {a.tmb_mut_per_mb:>8.2f} / "
        f"{b.tmb_mut_per_mb:>8.2f}  = {_safe_ratio(a.tmb_mut_per_mb, b.tmb_mut_per_mb)}"
    )
