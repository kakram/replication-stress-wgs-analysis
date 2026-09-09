"""
Canonical bibliography merge tool.

Per the invariants in CLAUDE.md (Bibliography Management):
  1. Westminster-Harvard parenthetical-year format only.
  2. Entries kept in alphabetical order by first-author surname, then by year.
  3. ALL merges from referencesDelta.docx must go through this script.
  4. references.docx is SHA-256-backed-up before any modification.

Usage:
    python scripts/merge_references.py [--references PATH] [--delta PATH] [--dry-run]
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import unicodedata
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

from docx import Document

REPO = Path(__file__).resolve().parent.parent
DEFAULT_REFS = REPO / "thesis" / "references.docx"
DEFAULT_DELTA = REPO / "thesis" / "referencesDelta.docx"
OUT = REPO / "outputs"
OUT.mkdir(exist_ok=True)


# ============================ shared parsing helpers ====================== #

def norm(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    s = s.replace(" ", " ").replace(" ", " ")
    return re.sub(r"\s+", " ", s).strip()


def strip_accents(s: str) -> str:
    nf = unicodedata.normalize("NFKD", s)
    return "".join(c for c in nf if not unicodedata.combining(c))


def make_key(surname: str, year: str) -> str:
    s = strip_accents(surname).lower()
    s = re.sub(r"[^a-z]", "", s)
    return f"{s}_{year}"


def sort_key_for(surname: str, year: str) -> tuple[str, str]:
    s = strip_accents(surname).lower()
    s = re.sub(r"[^a-z]", "", s)
    return (s, year)


# First character of a surname can be ASCII A-Z or an accented uppercase
# letter from Latin-1 Supplement / Latin Extended-A / Extended-B
# (covers Á, À, Â, Ã, Ä, Å, Ç, É, Í, Ñ, Ó, Ö, Ú, Ü, Ý, etc.).
# Note: \w in Python 3 already matches subsequent diacritic letters.
_LETTER_UPPER = r"A-ZÀ-ɏ"
_SURNAME = (
    rf"[{_LETTER_UPPER}][\w'’\-]+"
    rf"(?:\s+(?:van|de|der|den|von|la|le|du|del|di|da|dos|of)\s+"
    rf"[{_LETTER_UPPER}][\w'’\-]+)?"
)
SURNAME_ANCHORED = re.compile(rf"^({_SURNAME})")
YEAR_PARENS = re.compile(r"\(((?:19|20)\d{2}[a-z]?)\)")
# bare-year forms at the head: ", YYYY[, Month]." after the author block.
YEAR_BARE_HEAD = re.compile(
    r"^(?P<authors>.+?),\s*(?P<year>(?:19|20)\d{2}[a-z]?)"
    r"(?:,\s*[A-Z][a-z]+)?\s*\.\s*(?P<rest>.*)$"
)


# ============================ parsed entries ============================== #

@dataclass
class RefEntry:
    surname: str
    year: str
    key: str
    sort_key: tuple[str, str]
    paragraph: object  # docx.text.paragraph.Paragraph
    paragraph_index: int  # 0-based body paragraph index at parse time

    @property
    def text(self) -> str:
        return norm(self.paragraph.text)


@dataclass
class DeltaEntry:
    surname: Optional[str]
    year: Optional[str]
    key: Optional[str]
    sort_key: Optional[tuple[str, str]]
    source_para_index: int
    original_text: str
    normalised_text: str
    year_format: str          # parens | bare_year | bare_year_month | other
    parseable: bool
    deviations: list[str] = field(default_factory=list)
    italics_lost: bool = False


# ============================ parsing ===================================== #

def parse_ref_paragraph(text: str) -> Optional[tuple[str, str]]:
    """Permissive: extract (surname, year) from a Westminster-Harvard
    or near-Harvard reference entry. Returns None if not parseable."""
    t = norm(text)
    if not t or t.lower() == "references":
        return None
    m = SURNAME_ANCHORED.match(t)
    if not m:
        return None
    surname = m.group(1)
    y = YEAR_PARENS.search(t)
    if y:
        return surname, y.group(1)
    m_bare = YEAR_BARE_HEAD.match(t)
    if m_bare:
        return surname, m_bare.group("year")
    return None


def parse_references_doc(doc) -> tuple[list[RefEntry], list[int]]:
    """Return (parsed_entries, unparseable_paragraph_indexes)."""
    entries: list[RefEntry] = []
    unparseable: list[int] = []
    for i, p in enumerate(doc.paragraphs):
        text = (p.text or "").strip()
        if not text or text.lower() == "references":
            continue
        parsed = parse_ref_paragraph(p.text)
        if parsed is None:
            unparseable.append(i)
            continue
        surname, year = parsed
        entries.append(RefEntry(
            surname=surname,
            year=year,
            key=make_key(surname, year),
            sort_key=sort_key_for(surname, year),
            paragraph=p,
            paragraph_index=i,
        ))
    return entries, unparseable


# ============================ delta normalisation ========================= #

def detect_year_format(text: str) -> str:
    """Return one of: 'parens', 'bare_year_month', 'bare_year', 'other'."""
    t = norm(text)
    m = SURNAME_ANCHORED.match(t)
    if not m:
        return "other"
    # parens form: "Author (YYYY)" or "Author, X.Y. (YYYY)" - the year-parens
    # appears soon after the author block, not buried mid-sentence.
    # Heuristic: the FIRST (YYYY) found is within the first 250 chars AND
    # is preceded only by author-name characters (letters, commas, dots,
    # apostrophes, hyphens, ampersands, spaces, "and").
    yp = YEAR_PARENS.search(t)
    if yp and yp.start() < 250:
        head = t[:yp.start()]
        if re.fullmatch(r"[A-Za-zÀ-ÿ'’\-\.,&\s]+", head):
            return "parens"
    # bare-year-with-month: ", YYYY, Month."
    if re.match(r"^.+?,\s*(?:19|20)\d{2}[a-z]?,\s*[A-Z][a-z]+\s*\.", t):
        return "bare_year_month"
    # bare-year-with-period: ", YYYY."
    if re.match(r"^.+?,\s*(?:19|20)\d{2}[a-z]?\s*\.", t):
        return "bare_year"
    return "other"


def normalise_to_parens_year(text: str, year: str, year_format: str) -> str:
    """Rewrite text from bare-year form to Westminster-Harvard parens form."""
    t = norm(text)
    if year_format == "bare_year_month":
        return re.sub(
            rf",\s*{re.escape(year)},\s*[A-Z][a-z]+\s*\.",
            f" ({year}).",
            t, count=1,
        )
    if year_format == "bare_year":
        return re.sub(
            rf",\s*{re.escape(year)}\s*\.",
            f" ({year}).",
            t, count=1,
        )
    return t


def detect_style_deviations(text: str) -> list[str]:
    """Flag Westminster-Harvard style deviations the script will NOT auto-fix."""
    flags: list[str] = []
    if re.search(r"\w\s+&\s+\w", text):
        flags.append("uses '&' between authors (Westminster prefers 'and')")
    if '"' in text or "”" in text or "“" in text:
        flags.append("uses double quotes (Westminster prefers single quotes)")
    if re.search(r"\bno\.\s*\d", text):
        flags.append("uses 'no.' for issue number (MLA-style remnant)")
    return flags


def has_italics(p) -> bool:
    return any(run.italic for run in p.runs if run.italic)


def parse_delta_doc(doc) -> list[DeltaEntry]:
    """Parse all non-empty delta paragraphs, including unparseable ones."""
    out: list[DeltaEntry] = []
    for i, p in enumerate(doc.paragraphs):
        text = norm(p.text)
        if not text or text.lower() == "references":
            continue
        year_format = detect_year_format(text)
        parsed = parse_ref_paragraph(text)
        if parsed is None:
            out.append(DeltaEntry(
                surname=None, year=None, key=None, sort_key=None,
                source_para_index=i,
                original_text=text,
                normalised_text=text,
                year_format=year_format,
                parseable=False,
                deviations=detect_style_deviations(text),
                italics_lost=False,
            ))
            continue
        surname, year = parsed
        if year_format in ("bare_year", "bare_year_month"):
            new_text = normalise_to_parens_year(text, year, year_format)
            italics_lost = has_italics(p)
        else:
            new_text = text
            italics_lost = False
        out.append(DeltaEntry(
            surname=surname,
            year=year,
            key=make_key(surname, year),
            sort_key=sort_key_for(surname, year),
            source_para_index=i,
            original_text=text,
            normalised_text=new_text,
            year_format=year_format,
            parseable=True,
            deviations=detect_style_deviations(new_text),
            italics_lost=italics_lost,
        ))
    return out


# ============================ insertion =================================== #

def reference_paragraph_style(doc):
    for p in doc.paragraphs:
        t = (p.text or "").strip()
        if t and t.lower() != "references":
            return p.style
    return None


def find_insertion_target(new_sort_key: tuple[str, str],
                          ref_entries_sorted: list[RefEntry]) -> Optional[RefEntry]:
    """Return the first existing RefEntry whose sort_key > new_sort_key,
    or None if the new entry sorts after every existing entry."""
    for ref in ref_entries_sorted:
        if ref.sort_key > new_sort_key:
            return ref
    return None


def _copy_runs(source_paragraph, target_paragraph) -> None:
    for run in source_paragraph.runs:
        nr = target_paragraph.add_run(run.text)
        if run.italic is not None:
            nr.italic = run.italic
        if run.bold is not None:
            nr.bold = run.bold
        if run.underline is not None:
            nr.underline = run.underline


def insert_entry(doc, entry: DeltaEntry, source_para,
                 target_ref: Optional[RefEntry], ref_style) -> int:
    """Insert a new paragraph for `entry` either before `target_ref.paragraph`
    or at the end of the document. Returns the 0-based body paragraph index
    of the newly-inserted paragraph."""
    if target_ref is not None:
        new_p = target_ref.paragraph.insert_paragraph_before(style=ref_style)
    else:
        new_p = doc.add_paragraph()
        if ref_style is not None:
            new_p.style = ref_style

    if entry.parseable and entry.year_format == "parens" and not entry.italics_lost:
        # Conformant entry: preserve run-level formatting from source.
        _copy_runs(source_para, new_p)
    else:
        # Normalised or unparseable: emit as plain text (no italics).
        new_p.add_run(entry.normalised_text)

    # Look up the new paragraph's body index by element identity.
    elem = new_p._element
    for idx, p in enumerate(doc.paragraphs):
        if p._element is elem:
            return idx
    return -1


# ============================ hashing ===================================== #

def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# ============================ main pipeline =============================== #

def run(references_path: Path, delta_path: Path, dry_run: bool,
        report_dir: Path = OUT) -> int:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log: list[str] = []

    def say(msg: str) -> None:
        print(msg)
        log.append(msg)

    say("=" * 72)
    say(f"merge_references.py — {datetime.now().isoformat(timespec='seconds')}")
    say(f"  references : {references_path}")
    say(f"  delta      : {delta_path}")
    say(f"  dry-run    : {dry_run}")
    say("=" * 72)

    # ---- pre-flight ----
    if not references_path.exists():
        say(f"FATAL: references file not found at {references_path}")
        return 2
    if not delta_path.exists():
        say(f"FATAL: delta file not found at {delta_path}")
        return 2

    refs_doc = Document(str(references_path))
    delta_doc = Document(str(delta_path))
    n_ref_paras = len(refs_doc.paragraphs)
    n_delta_paras = len(delta_doc.paragraphs)
    n_ref_non_empty = sum(1 for p in refs_doc.paragraphs if p.text.strip())
    n_delta_non_empty = sum(1 for p in delta_doc.paragraphs if p.text.strip())
    say(f"references.docx     : {n_ref_paras} paragraphs "
        f"({n_ref_non_empty} non-empty)")
    say(f"referencesDelta.docx: {n_delta_paras} paragraphs "
        f"({n_delta_non_empty} non-empty)")

    # ---- backup + SHA-256 verify ----
    backup_path = references_path.with_name(
        f"{references_path.stem}_backup_{timestamp}{references_path.suffix}"
    )
    if dry_run:
        say(f"[dry-run] Backup would be written to {backup_path.name}")
    else:
        shutil.copy2(references_path, backup_path)
        orig_hash = sha256_of(references_path)
        backup_hash = sha256_of(backup_path)
        if orig_hash != backup_hash:
            say("FATAL: backup SHA-256 verification failed. Aborting.")
            say(f"  orig   : {orig_hash}")
            say(f"  backup : {backup_hash}")
            return 3
        say(f"Backup written: {backup_path.name}")
        say(f"SHA-256        : {orig_hash} (verified)")

    # ---- parse both ----
    ref_entries, ref_unparseable = parse_references_doc(refs_doc)
    delta_entries = parse_delta_doc(delta_doc)
    delta_parseable = [d for d in delta_entries if d.parseable]
    delta_unparseable = [d for d in delta_entries if not d.parseable]

    ref_keys = {r.key for r in ref_entries}
    delta_keys = {d.key for d in delta_parseable}
    shared_keys = ref_keys & delta_keys
    new_keys_set = delta_keys - ref_keys
    new_entries = [d for d in delta_parseable if d.key in new_keys_set]
    # Sort by alphabetical sort_key for deterministic insertion order.
    new_entries.sort(key=lambda d: d.sort_key)

    say(f"Parsed references   : {len(ref_entries)} entries "
        f"(unparseable paragraphs: {len(ref_unparseable)})")
    say(f"Parsed delta        : {len(delta_parseable)} entries "
        f"(unparseable paragraphs: {len(delta_unparseable)})")
    say(f"Shared keys         : {len(shared_keys)}")
    say(f"NEW keys to insert  : {len(new_entries)}")

    if delta_unparseable:
        say("Unparseable delta entries:")
        for d in delta_unparseable:
            say(f"  para {d.source_para_index}: {d.original_text[:120]}…")

    # ---- sort existing refs for insertion lookup ----
    ref_sorted = sorted(ref_entries, key=lambda r: r.sort_key)
    ref_style = reference_paragraph_style(refs_doc)

    # ---- plan insertions ----
    insertion_plan: list[dict] = []
    for entry in new_entries:
        target_ref = find_insertion_target(entry.sort_key, ref_sorted)
        if target_ref is None:
            location_desc = "end of document"
            target_idx = None
        else:
            target_idx = target_ref.paragraph_index
            location_desc = (f"before para {target_idx} "
                             f"({target_ref.surname} {target_ref.year})")
        insertion_plan.append({
            "entry": entry,
            "target_ref": target_ref,
            "location_desc": location_desc,
            "target_idx": target_idx,
        })

    # ---- apply (or skip in dry-run) ----
    inserted_at: dict[str, int] = {}
    if not dry_run:
        for step in insertion_plan:
            entry = step["entry"]
            src_para = delta_doc.paragraphs[entry.source_para_index]
            new_idx = insert_entry(refs_doc, entry, src_para,
                                   step["target_ref"], ref_style)
            inserted_at[entry.key] = new_idx
        refs_doc.save(str(references_path))
        # re-load to report after-state
        after_doc = Document(str(references_path))
        after_entries, _ = parse_references_doc(after_doc)
        n_after = len(after_entries)
        say(f"Wrote merged references.docx with {len(new_entries)} new entries.")
        say(f"After merge         : {n_after} entries")
    else:
        say("[dry-run] No file changes written.")
        n_after = len(ref_entries) + len(new_entries)
        say(f"After merge (proj.) : {n_after} entries")

    # ---- write report ----
    report_path = report_dir / f"merge_references_report_{timestamp}.md"
    write_report(
        report_path=report_path,
        timestamp=timestamp,
        dry_run=dry_run,
        references_path=references_path,
        delta_path=delta_path,
        backup_path=backup_path if not dry_run else None,
        n_ref_before=len(ref_entries),
        n_after=n_after,
        shared_keys=shared_keys,
        new_entries=new_entries,
        insertion_plan=insertion_plan,
        inserted_at=inserted_at,
        delta_unparseable=delta_unparseable,
        ref_unparseable=ref_unparseable,
        run_log=log,
    )
    try:
        report_display = report_path.relative_to(REPO)
    except ValueError:
        report_display = report_path
    say(f"Report: {report_display}")
    return 0


# ============================ report ====================================== #

def write_report(*, report_path: Path, timestamp: str, dry_run: bool,
                 references_path: Path, delta_path: Path,
                 backup_path: Optional[Path],
                 n_ref_before: int, n_after: int,
                 shared_keys: set,
                 new_entries: list[DeltaEntry],
                 insertion_plan: list[dict],
                 inserted_at: dict,
                 delta_unparseable: list[DeltaEntry],
                 ref_unparseable: list[int],
                 run_log: list[str]) -> None:
    with report_path.open("w", encoding="utf-8") as fh:
        fh.write("# Bibliography Merge Report\n\n")
        fh.write(f"- **Timestamp:** {timestamp}\n")
        fh.write(f"- **Mode:** {'DRY-RUN (no files written)' if dry_run else 'APPLY'}\n")
        fh.write(f"- **references.docx:** `{references_path}`\n")
        fh.write(f"- **referencesDelta.docx:** `{delta_path}`\n")
        if backup_path is not None:
            fh.write(f"- **Backup:** `{backup_path.name}` (SHA-256 verified)\n")
        fh.write(f"- **Entries before:** {n_ref_before}\n")
        fh.write(f"- **Entries after{' (projected)' if dry_run else ''}:** "
                 f"{n_after}\n")
        fh.write(f"- **Shared keys (skipped):** {len(shared_keys)}\n")
        fh.write(f"- **New entries inserted:** {len(new_entries)}\n\n")

        fh.write("## New entries — insertion plan\n\n")
        if not new_entries:
            fh.write("_No new entries to insert._\n\n")
        else:
            fh.write("| # | Surname | Year | Original delta format | "
                     "Normalised form | Insertion point | "
                     "Italics preserved? | Deviations flagged |\n")
            fh.write("|---|---------|------|-----------------------|"
                     "-----------------|-----------------|"
                     "--------------------|-------------------|\n")
            for i, step in enumerate(insertion_plan, 1):
                e = step["entry"]
                loc = step["location_desc"]
                if not dry_run and e.key in inserted_at:
                    loc = f"{loc} → wrote at body para {inserted_at[e.key]}"
                orig = e.original_text.replace("|", "\\|")
                if len(orig) > 110:
                    orig = orig[:107] + "..."
                norm_txt = e.normalised_text.replace("|", "\\|")
                if len(norm_txt) > 110:
                    norm_txt = norm_txt[:107] + "..."
                italics = ("no — plain text (italics lost)" if e.italics_lost
                           else "n/a (no italics in source)"
                           if e.year_format in ("bare_year", "bare_year_month")
                           else "yes (runs copied)")
                devs = "; ".join(e.deviations) if e.deviations else "—"
                fh.write(f"| {i} | {e.surname} | {e.year} | "
                         f"`{e.year_format}` — {orig} | "
                         f"{norm_txt} | {loc} | {italics} | {devs} |\n")
            fh.write("\n")

        fh.write("## Entries flagged for manual style review\n\n")
        flagged = [step["entry"] for step in insertion_plan
                   if step["entry"].deviations or step["entry"].italics_lost
                   or step["entry"].year_format == "other"]
        if not flagged:
            fh.write("_None._\n\n")
        else:
            fh.write("| Surname | Year | Reason |\n")
            fh.write("|---------|------|--------|\n")
            for e in flagged:
                reasons: list[str] = []
                if e.year_format == "other":
                    reasons.append("year-format not auto-normalisable "
                                   "(needs manual rewrite to Harvard)")
                if e.italics_lost:
                    reasons.append("italics lost during year-format "
                                   "normalisation (re-apply if needed)")
                reasons.extend(e.deviations)
                fh.write(f"| {e.surname or '—'} | {e.year or '—'} | "
                         f"{'; '.join(reasons)} |\n")
            fh.write("\n")

        fh.write("## Entries that could not be parsed\n\n")
        if not delta_unparseable:
            fh.write("_None._\n\n")
        else:
            fh.write("| Source para | Text |\n")
            fh.write("|------------:|------|\n")
            for d in delta_unparseable:
                t = d.original_text.replace("|", "\\|")
                if len(t) > 200:
                    t = t[:197] + "..."
                fh.write(f"| {d.source_para_index} | {t} |\n")
            fh.write("\n")

        if ref_unparseable:
            fh.write("## references.docx paragraphs the parser skipped\n\n")
            fh.write(f"_{len(ref_unparseable)} non-empty paragraph(s) in "
                     "`references.docx` could not be parsed as a reference. "
                     "These are NOT modified by this script._\n\n")

        fh.write("## Run log\n\n```\n")
        fh.write("\n".join(run_log))
        fh.write("\n```\n")


# ============================ CLI ========================================= #

def parse_cli(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--references", type=Path, default=DEFAULT_REFS,
                   help=f"references.docx (default: {DEFAULT_REFS})")
    p.add_argument("--delta", type=Path, default=DEFAULT_DELTA,
                   help=f"referencesDelta.docx (default: {DEFAULT_DELTA})")
    p.add_argument("--dry-run", action="store_true",
                   help="Analyse and write report only; do not modify "
                        "references.docx.")
    p.add_argument("--report-dir", type=Path, default=OUT,
                   help="Directory for the merge report markdown.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_cli(argv)
    return run(
        references_path=args.references,
        delta_path=args.delta,
        dry_run=args.dry_run,
        report_dir=args.report_dir,
    )


if __name__ == "__main__":
    sys.exit(main())
