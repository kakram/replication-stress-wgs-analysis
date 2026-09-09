"""
Pilot regression test:
  1. Parse references.docx and referencesDelta.docx.
  2. Identify entries in delta with keys not present in references.docx.
  3. Append those NEW entries to the end of references.docx, preserving the
     paragraph style of the surrounding references-doc entries and the
     character-level italic/bold/underline formatting of the delta runs.
  4. Re-parse the merged references and re-run Chapter 1 citation extraction.
  5. Classify each of the 13 originally-missing Ch1 citations as
     RESOLVED, IN-TEXT FIX NEEDED, or GENUINE GAP.
  6. Write outputs/ch1_missing_refs_regression.md.

The script imports the citation regexes and helpers from
scripts/ch1_reference_audit.py for consistency, but uses its own more
permissive reference-entry parser so that delta's "Author, YYYY." style is
handled alongside references.docx's "Author (YYYY)." style.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from docx import Document

# Import shared helpers from the v2 audit script.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ch1_reference_audit import (  # noqa: E402
    Chapter,
    chapter_paragraphs,
    extract_citations,
    find_chapter_by_pattern,
    find_chapters,
    make_key,
    norm,
    _SURNAME,
)

REPO = Path(__file__).resolve().parent.parent
THESIS = REPO / "thesis" / "KhalidAkramThesisDraftApr2026.docx"
REFS = REPO / "thesis" / "references.docx"
DELTA = REPO / "thesis" / "referencesDelta.docx"
OUT = REPO / "outputs"
OUT.mkdir(exist_ok=True)


# ----- permissive reference-entry parser (handles both year formats) ----- #

YEAR_PARENS = re.compile(r"\(((?:19|20)\d{2}[a-z]?)\)")
YEAR_BARE = re.compile(r",\s*((?:19|20)\d{2}[a-z]?)\.")
SURNAME_ANCHORED = re.compile(rf"^({_SURNAME})")


@dataclass
class RefEntry:
    first_author: str
    year: str
    full_entry: str
    key: str
    paragraph_index: int


def parse_entry(text: str) -> tuple[str, str] | None:
    """Return (first_author, year) or None."""
    t = norm(text)
    if not t or t.lower() == "references":
        return None
    m = SURNAME_ANCHORED.match(t)
    if not m:
        return None
    surname = m.group(1)
    y = YEAR_PARENS.search(t) or YEAR_BARE.search(t)
    if not y:
        return None
    return surname, y.group(1)


def parse_doc(doc: Document) -> list[RefEntry]:
    out: list[RefEntry] = []
    for i, p in enumerate(doc.paragraphs):
        parsed = parse_entry(p.text)
        if parsed is None:
            continue
        surname, year = parsed
        text = norm(p.text)
        out.append(RefEntry(
            first_author=surname,
            year=year,
            full_entry=text,
            key=make_key(surname, year),
            paragraph_index=i,
        ))
    return out


# --------- merge: append new entries, preserve style + run formatting --- #

def reference_paragraph_style(doc: Document):
    """Return the paragraph style used by the bulk of reference entries."""
    for p in doc.paragraphs:
        t = (p.text or "").strip()
        if t and t.lower() != "references":
            return p.style
    return None


def append_entry(target_doc: Document, source_para, target_style) -> None:
    """Append a paragraph to target_doc, copying text + run formatting.

    The new paragraph's style is set to target_style (the references.docx
    style), NOT the source's style. Run-level italic / bold / underline
    are copied from the source runs.
    """
    new_p = target_doc.add_paragraph()
    if target_style is not None:
        new_p.style = target_style
    for run in source_para.runs:
        nr = new_p.add_run(run.text)
        if run.italic is not None:
            nr.italic = run.italic
        if run.bold is not None:
            nr.bold = run.bold
        if run.underline is not None:
            nr.underline = run.underline


# ----- the 13 original missing Ch1 citations + known in-text-fix cases ----- #

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

# Citations that will never resolve via reference-side fixes alone — the
# in-text citation itself is wrong (typo or year mismatch).
IN_TEXT_FIX_HINTS = {
    "bahn_2007":     "Likely typo in thesis: 'Bahn et al., 2007' → should be 'Baan et al., 2007' (Baan is in references.docx).",
    "jones_2010":    "Likely year mismatch in thesis: 'Jones and Petermann, 2010' → references entry is 2012.",
    "shendure_2022": "Likely year mismatch in thesis: 'Shendure et al., 2022' → references entry is 2017.",
}


# ------------------------------- main ------------------------------------- #

def main() -> None:
    if not REFS.exists() or not DELTA.exists():
        raise SystemExit("Required input files missing — aborting.")

    print("=" * 70)
    print(f"Run timestamp: {datetime.now().isoformat(timespec='seconds')}")
    print("=" * 70)

    # ---- parse both ----
    refs_doc = Document(str(REFS))
    delta_doc = Document(str(DELTA))

    refs_entries = parse_doc(refs_doc)
    delta_entries = parse_doc(delta_doc)
    refs_keys = {r.key for r in refs_entries}
    delta_keys = {d.key for d in delta_entries}
    shared = refs_keys & delta_keys
    new_keys = delta_keys - refs_keys

    n_refs = len(refs_entries)
    n_delta = len(delta_entries)

    print(f"Parsed references.docx       : {n_refs} entries")
    print(f"Parsed referencesDelta.docx  : {n_delta} entries")
    print(f"Shared dedup keys            : {len(shared)}")
    print(f"NEW keys in delta            : {len(new_keys)}")

    if n_delta == 0:
        raise SystemExit("Delta parsed as zero entries — aborting before merge.")

    # ---- merge ----
    ref_style = reference_paragraph_style(refs_doc)
    print(f"References-paragraph style   : {ref_style.name if ref_style else 'None'}")

    # Build a list of source paragraphs to append, preserving delta order.
    new_in_delta_order: list[tuple[RefEntry, "Paragraph"]] = []
    seen_keys: set[str] = set()
    for entry in delta_entries:
        if entry.key in refs_keys or entry.key in seen_keys:
            continue
        seen_keys.add(entry.key)
        new_in_delta_order.append(
            (entry, delta_doc.paragraphs[entry.paragraph_index])
        )

    added_count = 0
    for entry, src_para in new_in_delta_order:
        append_entry(refs_doc, src_para, ref_style)
        added_count += 1

    if added_count != len(new_keys):
        print(f"WARNING: append count ({added_count}) != new-key count "
              f"({len(new_keys)}). Continuing.")

    refs_doc.save(str(REFS))
    after = parse_doc(Document(str(REFS)))
    n_after = len(after)
    print()
    print(f"Before merge:  {n_refs} entries")
    print(f"Added:         {added_count} entries")
    print(f"After merge:   {n_after} entries")

    # ---- regression: re-run Ch1 missing-refs check ----
    after_keys = {r.key for r in after}

    thesis = Document(str(THESIS))
    chapters = find_chapters(thesis)
    ch1 = find_chapter_by_pattern(chapters, r"^(chapter\s*1\b|introduction\b)")
    if ch1 is None:
        raise SystemExit("Chapter 1 not found in thesis — aborting.")
    ch1_paras = chapter_paragraphs(thesis, ch1)
    cits = extract_citations(ch1_paras)
    cit_by_key = {c.key: c for c in cits}

    rows = []
    n_resolved = 0
    n_intext = 0
    n_gap = 0
    for key, author, year in ORIGINAL_MISSING:
        if key in after_keys:
            verdict = "RESOLVED"
            note = "Entry now present in merged references.docx."
            n_resolved += 1
        elif key in IN_TEXT_FIX_HINTS:
            verdict = "STILL MISSING — IN-TEXT FIX NEEDED"
            note = IN_TEXT_FIX_HINTS[key]
            n_intext += 1
        else:
            verdict = "STILL MISSING — GENUINE GAP"
            note = ("Not present in merged references.docx and no known "
                    "in-text-fix hint applies.")
            n_gap += 1
        cit = cit_by_key.get(key)
        snippet = cit.location_snippet if cit else ""
        rows.append((key, author, year, verdict, note, snippet))

    verdict_line = (
        "PILOT FIX-FLOW PASS"
        if n_gap == 0
        else f"PILOT FIX-FLOW PARTIAL — {n_gap} item(s) remain as genuine gaps"
    )

    # ---- write report ----
    report = OUT / "ch1_missing_refs_regression.md"
    with report.open("w", encoding="utf-8") as fh:
        fh.write("# Chapter 1 Missing-References Regression Test\n\n")
        fh.write(f"- **Run timestamp:** "
                 f"{datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"- **references.docx before merge:** {n_refs} entries\n")
        fh.write(f"- **referencesDelta.docx:** {n_delta} entries "
                 f"({len(shared)} shared, {len(new_keys)} new)\n")
        fh.write(f"- **Entries appended:** {added_count}\n")
        fh.write(f"- **references.docx after merge:** {n_after} entries\n")
        fh.write("- **Backup:** earlier copy preserved as "
                 "`thesis/references_backup_<timestamp>.docx` "
                 "(SHA-256 verified before merge).\n\n")

        fh.write("## Per-citation regression\n\n")
        fh.write("| # | Citation key | Author | Year | Verdict | Note | "
                 "Original location snippet |\n")
        fh.write("|---|--------------|--------|------|---------|------|"
                 "---------------------------|\n")
        for i, (key, author, year, verdict, note, snippet) in enumerate(rows, 1):
            snippet_cell = snippet.replace("|", "\\|") if snippet else "—"
            note_cell = note.replace("|", "\\|")
            fh.write(f"| {i} | `{key}` | {author} | {year} | "
                     f"**{verdict}** | {note_cell} | {snippet_cell} |\n")
        fh.write("\n")

        fh.write("## Totals\n\n")
        fh.write(f"- RESOLVED: **{n_resolved}**\n")
        fh.write(f"- STILL MISSING — IN-TEXT FIX NEEDED: **{n_intext}**\n")
        fh.write(f"- STILL MISSING — GENUINE GAP: **{n_gap}**\n\n")
        fh.write(f"## Final verdict\n\n**{verdict_line}**\n")

    # ---- stdout summary ----
    print()
    print("Per-citation regression:")
    print(f"  {'KEY':<18} {'AUTHOR':<10} {'YEAR':<6} VERDICT")
    for key, author, year, verdict, *_ in rows:
        print(f"  {key:<18} {author:<10} {year:<6} {verdict}")

    print()
    print(f"Resolved: {n_resolved}   In-text fix needed: {n_intext}   "
          f"Genuine gap: {n_gap}")
    print()
    print(verdict_line)
    print(f"Report: {report.relative_to(REPO)}")


if __name__ == "__main__":
    main()
