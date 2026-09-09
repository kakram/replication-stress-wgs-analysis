"""
Narrow Chapter-1 missing-references regression check.

Re-runs Task 3 Section A (missing references) of the audit, scoped to
Chapter 1 only, against the current state of references.docx. Does NOT
merge or modify any files.

Output: outputs/ch1_missing_refs_regression_final.md
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

from docx import Document

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ch1_reference_audit import (  # noqa: E402
    chapter_paragraphs,
    extract_citations,
    find_chapter_by_pattern,
    find_chapters,
)
from merge_references import parse_references_doc  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
THESIS = REPO / "thesis" / "KhalidAkramThesisDraftApr2026.docx"
REFS = REPO / "thesis" / "references.docx"
OUT = REPO / "outputs"


ORIGINAL_MISSING = [
    ("bahn_2007",     "Bahn",     "2007"),
    ("banerji_2012",  "Banerji",  "2012"),
    ("bucknall_1973", "Bucknall", "1973"),
    ("diffley_2011",  "Diffley",  "2011"),
    ("ellis_2013",    "Ellis",    "2013"),
    ("grissa_2007",   "Grissa",   "2007"),
    ("grissa_2008",   "Grissa",   "2008"),
    ("iliakis_1982",  "Iliakis",  "1982"),
    ("ishino_1987",   "Ishino",   "1987"),
    ("jones_2010",    "Jones",    "2010"),
    ("oberle_1991",   "Oberle",   "1991"),
    ("shendure_2022", "Shendure", "2022"),
    ("solaiman_2021", "Solaiman", "2021"),
]

IN_TEXT_FIX_HINTS = {
    "bahn_2007":     "Likely typo in thesis: 'Bahn et al., 2007' → should be 'Baan et al., 2007' (Baan is in references.docx).",
    "jones_2010":    "Likely year mismatch in thesis: 'Jones and Petermann, 2010' → references entry is 2012.",
    "shendure_2022": "Likely year mismatch in thesis: 'Shendure et al., 2022' → references entry is 2017.",
}


def main() -> int:
    refs_doc = Document(str(REFS))
    thesis = Document(str(THESIS))

    ref_entries, _ = parse_references_doc(refs_doc)
    ref_keys = {r.key for r in ref_entries}

    chapters = find_chapters(thesis)
    ch1 = find_chapter_by_pattern(chapters, r"^(chapter\s*1\b|introduction\b)")
    if ch1 is None:
        raise SystemExit("Chapter 1 heading not found.")
    cits = extract_citations(chapter_paragraphs(thesis, ch1))
    cit_by_key = {c.key: c for c in cits}

    rows = []
    n_res = n_intext = n_gap = 0
    for key, author, year in ORIGINAL_MISSING:
        if key in ref_keys:
            verdict = "RESOLVED"
            note = "Entry now present in references.docx."
            n_res += 1
        elif key in IN_TEXT_FIX_HINTS:
            verdict = "STILL MISSING — IN-TEXT FIX NEEDED"
            note = IN_TEXT_FIX_HINTS[key]
            n_intext += 1
        else:
            verdict = "STILL MISSING — GENUINE GAP"
            note = ("Not present in references.docx and no known "
                    "in-text-fix hint applies.")
            n_gap += 1
        snippet = cit_by_key.get(key)
        snippet_text = snippet.location_snippet if snippet else ""
        rows.append((key, author, year, verdict, note, snippet_text))

    verdict_line = (
        "PILOT FIX-FLOW PASS"
        if n_gap == 0
        else f"PILOT FIX-FLOW PARTIAL — {n_gap} item(s) remain as genuine gaps"
    )

    report = OUT / "ch1_missing_refs_regression_final.md"
    with report.open("w", encoding="utf-8") as fh:
        fh.write("# Chapter 1 Missing-References Regression — Final\n\n")
        fh.write(f"- **Run timestamp:** "
                 f"{datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"- **references.docx parsed entries:** {len(ref_entries)} "
                 "(after diacritic-led-surname parser fix)\n")
        fh.write(f"- **Chapter 1 paragraphs scanned:** "
                 f"{ch1.start}-{ch1.end - 1}\n")
        fh.write(f"- **Unique Ch1 in-text citations:** {len(cits)}\n\n")

        fh.write("## Per-citation regression\n\n")
        fh.write("| # | Citation key | Author | Year | Verdict | Note | "
                 "Original location snippet |\n")
        fh.write("|---|--------------|--------|------|---------|------|"
                 "---------------------------|\n")
        for i, (key, author, year, verdict, note, snippet) in enumerate(rows, 1):
            snip = snippet.replace("|", "\\|") if snippet else "—"
            note_cell = note.replace("|", "\\|")
            fh.write(f"| {i} | `{key}` | {author} | {year} | "
                     f"**{verdict}** | {note_cell} | {snip} |\n")
        fh.write("\n")

        fh.write("## Totals\n\n")
        fh.write(f"- RESOLVED: **{n_res}**\n")
        fh.write(f"- STILL MISSING — IN-TEXT FIX NEEDED: **{n_intext}**\n")
        fh.write(f"- STILL MISSING — GENUINE GAP: **{n_gap}**\n\n")

        fh.write(f"## Final verdict\n\n**{verdict_line}**\n")

    # stdout
    print(f"references.docx entries: {len(ref_entries)}")
    print(f"Ch1 citations: {len(cits)}\n")
    print(f"{'KEY':<18} {'AUTHOR':<10} {'YEAR':<6} VERDICT")
    for key, author, year, verdict, *_ in rows:
        print(f"{key:<18} {author:<10} {year:<6} {verdict}")
    print()
    print(f"Resolved: {n_res}   In-text fix needed: {n_intext}   "
          f"Genuine gap: {n_gap}")
    print()
    print(verdict_line)
    print(f"Report: {report.relative_to(REPO)}")
    return 0 if n_gap == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
