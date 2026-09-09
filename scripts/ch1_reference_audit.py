"""
Chapter 1 reference-checking pilot (v2).

v2 changes vs v1:
- Task 4 false-positive filter tightened (figure-caption / list-bullet / short
  sentence / paragraphs styled "Caption" are excluded).
- Task 4 no longer mechanically auto-generates topic summaries; instead the
  script writes outputs/_ch1_surviving_claims.json for Claude to consume and
  hand-craft topic summaries into outputs/ch1_citation_opportunities.md.
- Task 5 added: whole-thesis sweep for residual non-Harvard citation
  patterns (numeric brackets, paren digits, superscript runs).
- Chapter-boundary detection factored out into find_chapters() /
  find_chapter_by_pattern() so chapters 2-8 can reuse it later.
"""

from __future__ import annotations

import csv
import json
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path

from docx import Document

REPO = Path(__file__).resolve().parent.parent
THESIS = REPO / "thesis" / "KhalidAkramThesisDraftApr2026.docx"
REFS = REPO / "thesis" / "references.docx"
OUT = REPO / "outputs"
OUT.mkdir(exist_ok=True)


# ----------------------------- text helpers -------------------------------- #

def norm(s: str) -> str:
    s = unicodedata.normalize("NFKC", s)
    s = s.replace(" ", " ").replace(" ", " ")
    return re.sub(r"\s+", " ", s).strip()


def strip_accents(s: str) -> str:
    nf = unicodedata.normalize("NFKD", s)
    return "".join(c for c in nf if not unicodedata.combining(c))


def make_key(first_author: str, year: str) -> str:
    surname = strip_accents(first_author).lower()
    surname = re.sub(r"[^a-z]", "", surname)
    return f"{surname}_{year}"


# ------------------------- chapter boundary helpers ------------------------ #

@dataclass
class Chapter:
    heading: str
    start: int  # paragraph index of the heading paragraph itself
    end: int    # exclusive end paragraph index


def find_chapters(doc: Document) -> list[Chapter]:
    """Return all Heading 1 spans as Chapter records.

    Reusable for chapters 2-8 — pass the resulting list to
    find_chapter_by_pattern() with the appropriate regex.
    """
    h1s: list[tuple[int, str]] = []
    for i, p in enumerate(doc.paragraphs):
        style = (p.style.name if p.style else "") or ""
        text = norm(p.text)
        if "Heading 1" in style and text:
            h1s.append((i, text))
    chapters: list[Chapter] = []
    for j, (idx, name) in enumerate(h1s):
        end = h1s[j + 1][0] if j + 1 < len(h1s) else len(doc.paragraphs)
        chapters.append(Chapter(heading=name, start=idx, end=end))
    return chapters


def find_chapter_by_pattern(chapters: list[Chapter], pattern: str) -> Chapter | None:
    pat = re.compile(pattern, re.I)
    for c in chapters:
        if pat.search(c.heading):
            return c
    return None


# ----------------------------- paragraph model ----------------------------- #

@dataclass
class Para:
    idx: int
    style: str
    text: str


def chapter_paragraphs(doc: Document, chapter: Chapter) -> list[Para]:
    out: list[Para] = []
    for i in range(chapter.start, chapter.end):
        p = doc.paragraphs[i]
        out.append(Para(idx=i, style=(p.style.name or ""), text=norm(p.text)))
    return out


# ----------------------------- citation regexes ---------------------------- #

# First character of a surname can be ASCII A-Z or an accented uppercase
# letter from Latin-1 Supplement / Latin Extended-A / Extended-B
# (covers Á, À, Â, Ã, Ä, Å, Ç, É, Í, Ñ, Ó, Ö, Ú, Ü, Ý, etc.).
_LETTER_UPPER = r"A-ZÀ-ɏ"
_SURNAME = (
    rf"[{_LETTER_UPPER}][\w'’\-]+"
    rf"(?:\s+(?:van|de|der|den|von|la|le|du|del|di|da|dos|of)\s+"
    rf"[{_LETTER_UPPER}][\w'’\-]+)?"
)
_AUTHOR_GROUP = (
    rf"{_SURNAME}"
    rf"(?:\s+(?:and|&)\s+{_SURNAME})?"
    rf"(?:\s+et\s+al\.?)?"
)

PAREN_CIT = re.compile(
    rf"\(\s*((?:{_AUTHOR_GROUP},?\s*\d{{4}}[a-z]?(?:\s*[;,]\s*)?)+)\s*\)"
)
SINGLE_PAREN = re.compile(
    rf"({_AUTHOR_GROUP}),?\s*(\d{{4}}[a-z]?)"
)
NARR_CIT_NAMED = re.compile(
    rf"(?<![A-Za-z])({_AUTHOR_GROUP})\s*\((\d{{4}}[a-z]?)\)"
)


def extract_first_author(author_part: str) -> str:
    a = author_part.strip()
    a = re.sub(r"\s+et\s+al\.?$", "", a, flags=re.I)
    a = re.split(r"\s+(?:and|&)\s+", a, maxsplit=1)[0]
    return a.strip().rstrip(",")


@dataclass
class Citation:
    key: str
    first_author: str
    year: str
    location_snippet: str
    para_idx: int


def extract_citations(paras: list[Para]) -> list[Citation]:
    seen: dict[str, Citation] = {}
    for p in paras:
        text = p.text
        if not text:
            continue
        for m in PAREN_CIT.finditer(text):
            inside = m.group(1)
            for sm in SINGLE_PAREN.finditer(inside):
                first = extract_first_author(sm.group(1))
                year = sm.group(2)
                key = make_key(first, year)
                if key in seen:
                    continue
                s0 = max(0, m.start() - 60)
                s1 = min(len(text), m.end() + 60)
                seen[key] = Citation(key, first, year, text[s0:s1].strip(), p.idx)

        for m in NARR_CIT_NAMED.finditer(text):
            first = extract_first_author(m.group(1))
            if first.lower() in {
                "in", "the", "an", "as", "for", "by", "from", "see",
                "fig", "figure", "table", "e.g", "i.e", "et", "of",
            }:
                continue
            year = m.group(2)
            key = make_key(first, year)
            if key in seen:
                continue
            s0 = max(0, m.start() - 60)
            s1 = min(len(text), m.end() + 60)
            seen[key] = Citation(key, first, year, text[s0:s1].strip(), p.idx)
    return sorted(seen.values(), key=lambda c: (c.first_author.lower(), c.year))


# ----------------------------- references parsing -------------------------- #

REF_HEAD = re.compile(rf"^({_SURNAME})\b[^()]*?\((\d{{4}}[a-z]?)\)")


@dataclass
class Reference:
    first_author: str
    year: str
    full_entry: str
    key: str


def parse_references(doc: Document) -> list[Reference]:
    refs: list[Reference] = []
    for p in doc.paragraphs:
        text = norm(p.text)
        if not text or text.lower() == "references":
            continue
        m = REF_HEAD.match(text)
        if m:
            fa, yr = m.group(1), m.group(2)
        else:
            m2 = re.match(rf"^({_SURNAME})", text)
            y = re.search(r"\((\d{4}[a-z]?)\)", text)
            if not (m2 and y):
                continue
            fa, yr = m2.group(1), y.group(1)
        refs.append(Reference(fa, yr, text, make_key(fa, yr)))
    return refs


# ------------------------------- subsections ------------------------------- #

@dataclass
class Subsection:
    heading: str
    paras: list[Para]


def split_subsections(paras: list[Para]) -> list[Subsection]:
    subs: list[Subsection] = []
    cur = Subsection("(chapter preamble)", [])
    for p in paras:
        if "Heading 2" in p.style and p.text:
            if cur.paras:
                subs.append(cur)
            cur = Subsection(p.text, [])
        elif "Heading 1" in p.style:
            continue
        else:
            cur.paras.append(p)
    if cur.paras:
        subs.append(cur)
    return subs


# ---------------------- Task 4: claim-detection filters -------------------- #

HAS_CITATION = re.compile(
    rf"\(\s*[{_LETTER_UPPER}][\w'’\-]+[^()]*\d{{4}}"
    rf"|\b[{_LETTER_UPPER}][\w'’\-]+\s+"
    rf"(?:et al\.?\s*|and\s+[{_LETTER_UPPER}][\w'’\-]+\s+)?\(\d{{4}}\)"
)
CLAIM_CUES = re.compile(
    r"\b(studies|research|shown|evidence|reported|demonstrated|"
    r"associated with|linked to|known to|increase[sd]?|decrease[sd]?|"
    r"reduce[sd]?|elevated|approximately|around|over\s+\d|recent|recently|"
    r"emerging|leading cause|major cause|first identified|discovered|"
    r"in \d{4}|since \d{4}|by \d{4}|following the|widely|commonly|"
    r"frequently|estimated|roughly|whereas)\b",
    re.I,
)
PERCENT_OR_NUM = re.compile(r"\b\d+(?:\.\d+)?\s*%|\b\d{2,}\b")

# v2 exclusions
EXCLUDE_PREFIX = re.compile(
    r"^(fig\b|figure\b|left\s*:|right\s*:|top\s*:|bottom\s*:|"
    r"a\s*:|b\s*:|c\s*:|d\s*:)",
    re.I,
)
EXCLUDE_SUBSTR = re.compile(
    r"shown in fig|shown in figure|adapted from|\(see figure|"
    r"\bsee fig\b|\bsee figure\b|below\)|above\)|schematic of|"
    r"illustrates|depicts",
    re.I,
)


def is_caption(p: Para) -> bool:
    s = p.style or ""
    return s == "Caption" or s.startswith("Caption")


def find_uncited_claims_v2(paras: list[Para], max_per_sub: int = 3) -> list[str]:
    out: list[str] = []
    for p in paras:
        if not p.text or is_caption(p):
            continue
        sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z])", p.text)
        for s in sentences:
            s = s.strip()
            if len(s.split()) < 8:
                continue
            if len(s) < 60 or len(s) > 350:
                continue
            if HAS_CITATION.search(s):
                continue
            if EXCLUDE_PREFIX.search(s):
                continue
            if EXCLUDE_SUBSTR.search(s):
                continue
            if not (CLAIM_CUES.search(s) or PERCENT_OR_NUM.search(s)):
                continue
            out.append(s)
            if len(out) >= max_per_sub:
                return out
    return out


def find_uncited_claims_v1(paras: list[Para], max_per_sub: int = 3) -> list[str]:
    """v1 behaviour, kept only so we can quantify the diff after refactor."""
    out: list[str] = []
    for p in paras:
        if not p.text:
            continue
        sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z])", p.text)
        for s in sentences:
            s = s.strip()
            if len(s) < 60 or len(s) > 350:
                continue
            if HAS_CITATION.search(s):
                continue
            if CLAIM_CUES.search(s) or PERCENT_OR_NUM.search(s):
                if re.match(r"^(figure|fig\.|table|chapter)\s", s, re.I):
                    continue
                out.append(s)
                if len(out) >= max_per_sub:
                    return out
    return out


# ------------------ Task 5: residual numbered-citation sweep --------------- #

RE_BRACKET = re.compile(r"\[\d+(?:\s*[,\-–]\s*\d+)*\]")
RE_PAREN_DIGIT = re.compile(r"(?<!\w)\(\d{1,3}\)(?!\d)")
PAREN_DIGIT_BLOCKER = re.compile(
    r"(?:figure|fig\.?|table|tab\.?|eq\.?|equation|section|chapter|step)\s*$",
    re.I,
)
UNICODE_SUPER = "⁰¹²³⁴⁵⁶⁷⁸⁹"
RE_UNICODE_SUPER = re.compile(f"[{UNICODE_SUPER}]+")


def context_snippet(text: str, start: int, end: int,
                    before: int = 50, after: int = 50) -> str:
    s = max(0, start - before)
    e = min(len(text), end + after)
    return text[s:e].replace("\n", " ").strip()


def chapter_for_para(idx: int, chapters: list[Chapter]) -> str:
    for c in chapters:
        if c.start <= idx < c.end:
            return c.heading
    return "(front matter / unattached)"


def scan_numbered_residues(doc: Document,
                           chapters: list[Chapter]) -> list[dict]:
    findings: list[dict] = []
    for i, p in enumerate(doc.paragraphs):
        text = norm(p.text)
        chap = chapter_for_para(i, chapters)

        if text:
            for m in RE_BRACKET.finditer(text):
                findings.append({
                    "chapter": chap,
                    "para": i,
                    "pattern": m.group(0),
                    "context": context_snippet(text, m.start(), m.end()),
                })
            for m in RE_PAREN_DIGIT.finditer(text):
                pre = text[max(0, m.start() - 15):m.start()].rstrip()
                if PAREN_DIGIT_BLOCKER.search(pre):
                    continue
                findings.append({
                    "chapter": chap,
                    "para": i,
                    "pattern": m.group(0),
                    "context": context_snippet(text, m.start(), m.end()),
                })
            for m in RE_UNICODE_SUPER.finditer(text):
                findings.append({
                    "chapter": chap,
                    "para": i,
                    "pattern": f"unicode-super:{m.group(0)}",
                    "context": context_snippet(text, m.start(), m.end()),
                })

        # run-level superscript-numeric — catches Word's superscript style.
        raw = p.text or ""
        for run in p.runs:
            rt = (run.text or "").strip()
            if not rt:
                continue
            if not run.font.superscript:
                continue
            if not re.fullmatch(r"\d+(?:[,\-–]\d+)*", rt):
                continue
            pos = raw.find(run.text)
            if pos < 0:
                ctx = raw[:120]
            else:
                ctx = context_snippet(raw, pos, pos + len(run.text))
            findings.append({
                "chapter": chap,
                "para": i,
                "pattern": f"superscript:{rt}",
                "context": ctx,
            })
    return findings


# ----------------------------- audit reports ------------------------------- #

def write_audit_report(path: Path, ch1: Chapter, cits: list[Citation],
                       refs: list[Reference],
                       missing: list[Citation],
                       unused: list[Reference],
                       density_rows: list[tuple]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("# Chapter 1 Reference Audit\n\n")
        fh.write(f"- Chapter 1 paragraphs scanned: {ch1.start}-{ch1.end - 1} "
                 f"({ch1.end - ch1.start} paragraphs)\n")
        fh.write(f"- Unique in-text citations found: {len(cits)}\n")
        fh.write(f"- References parsed: {len(refs)}\n\n")

        fh.write("## A. Missing References\n\n")
        fh.write("Citations present in Chapter 1 with no matching entry in "
                 "references.docx (matched on first-author surname + year, "
                 "case- and accent-insensitive).\n\n")
        if not missing:
            fh.write("_None — every Chapter 1 citation has a matching "
                     "reference._\n\n")
        else:
            fh.write("| # | Citation | Year | Location snippet |\n")
            fh.write("|---|----------|------|------------------|\n")
            for i, c in enumerate(missing, 1):
                snip = c.location_snippet.replace("|", "\\|")
                fh.write(f"| {i} | {c.first_author} | {c.year} | {snip} |\n")
            fh.write("\n")

        fh.write("## B. Unused References\n\n")
        fh.write("> **Caveat:** these are entries in `references.docx` not "
                 "cited in Chapter 1. Many will legitimately be cited in "
                 "later chapters (Methods, Chapters 4-8). Treat this list "
                 "as a Chapter-1-scope hint, not a deletion list.\n\n")
        if not unused:
            fh.write("_None._\n\n")
        else:
            fh.write(f"Total unused (in Ch1 scope): **{len(unused)}**\n\n")
            fh.write("| # | First author | Year | Full entry |\n")
            fh.write("|---|--------------|------|------------|\n")
            for i, r in enumerate(unused, 1):
                entry = r.full_entry.replace("|", "\\|")
                if len(entry) > 220:
                    entry = entry[:217] + "..."
                fh.write(f"| {i} | {r.first_author} | {r.year} | {entry} |\n")
            fh.write("\n")

        fh.write("## C. Citation Density Map\n\n")
        fh.write("Density expressed as citations per 500 words. Subsections "
                 "with <3 per 500 words **and** ≥200 words are flagged.\n\n")
        fh.write("| Subsection | Words | Citations | per 500 w | Flag |\n")
        fh.write("|------------|------:|----------:|----------:|------|\n")
        for heading, wc, n, per500, flag in density_rows:
            fh.write(f"| {heading} | {wc} | {n} | {per500:.2f} | {flag} |\n")
        fh.write("\n")


def write_numbered_residues(path: Path, residues: list[dict]) -> None:
    by_chap: dict[str, list[dict]] = {}
    for r in residues:
        by_chap.setdefault(r["chapter"], []).append(r)

    with path.open("w", encoding="utf-8") as fh:
        fh.write("# Numbered / Footnoted Citation Residues "
                 "— Whole-Thesis Sweep\n\n")
        fh.write("Global sweep for residual non-Harvard citation patterns: "
                 "numeric brackets (e.g. `[12]`, `[1-3]`, `[1,2,5]`), "
                 "standalone paren-digit citations `(12)` not preceded by "
                 "Figure / Table / Eq. / Section / Chapter / Step, Unicode "
                 "superscript-digit sequences, and Word runs flagged as "
                 "`font.superscript=True` with purely numeric text.\n\n")

        n_idx = 0
        for chap_name in sorted(by_chap.keys()):
            fh.write(f"## {chap_name}\n\n")
            fh.write(f"_{len(by_chap[chap_name])} finding(s)_\n\n")
            fh.write("| # | Paragraph | Pattern | Context |\n")
            fh.write("|---|----------:|---------|---------|\n")
            for item in by_chap[chap_name]:
                n_idx += 1
                ctx = item["context"].replace("|", "\\|")
                pat = item["pattern"].replace("|", "\\|")
                fh.write(f"| {n_idx} | {item['para']} | `{pat}` | …{ctx}… |\n")
            fh.write("\n")

        fh.write("---\n\n")
        fh.write(f"**Summary:** {len(residues)} residual numbered citation"
                 f"(s) found across {len(by_chap)} chapter(s).\n")


# --------------------------------- main ------------------------------------ #

def main() -> None:
    thesis = Document(str(THESIS))
    refs_doc = Document(str(REFS))

    chapters = find_chapters(thesis)
    ch1 = find_chapter_by_pattern(chapters, r"^(chapter\s*1\b|introduction\b)")
    if ch1 is None:
        raise RuntimeError("Chapter 1 heading not found.")
    paras = chapter_paragraphs(thesis, ch1)
    print(f"Chapter 1: paragraph {ch1.start}-{ch1.end - 1} "
          f"({ch1.end - ch1.start} paragraphs)")

    # ---- Task 1 ----
    cits = extract_citations(paras)
    cits_csv = OUT / "ch1_citations_in_text.csv"
    with cits_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["citation_key", "first_author", "year", "location_snippet"])
        for c in cits:
            w.writerow([c.key, c.first_author, c.year, c.location_snippet])
    print(f"Task 1: {len(cits)} unique citations -> {cits_csv.name}")

    # ---- Task 2 ----
    refs = parse_references(refs_doc)
    refs_csv = OUT / "ch1_references_parsed.csv"
    with refs_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["first_author", "year", "full_entry"])
        for r in refs:
            w.writerow([r.first_author, r.year, r.full_entry])
    print(f"Task 2: {len(refs)} references -> {refs_csv.name}")

    # ---- Task 3 ----
    ref_keys = {r.key for r in refs}
    cit_keys = {c.key for c in cits}
    missing = [c for c in cits if c.key not in ref_keys]
    unused = [r for r in refs if r.key not in cit_keys]

    subs = split_subsections(paras)
    density_rows: list[tuple] = []
    for s in subs:
        wc = sum(len(p.text.split()) for p in s.paras)
        scs = extract_citations(s.paras)
        per500 = (len(scs) / wc * 500) if wc else 0.0
        flag = "POTENTIALLY UNDER-CITED" if (wc >= 200 and per500 < 3) else ""
        density_rows.append((s.heading, wc, len(scs), per500, flag))

    audit_md = OUT / "ch1_reference_audit.md"
    write_audit_report(audit_md, ch1, cits, refs, missing, unused, density_rows)
    print(f"Task 3: missing={len(missing)} unused={len(unused)} "
          f"flagged_subs={sum(1 for r in density_rows if r[4])} -> {audit_md.name}")

    # ---- Task 4 (v2) ----
    sub_to_claims: list[dict] = []
    v1_total = 0
    v2_total = 0
    for s in subs:
        v1_claims = find_uncited_claims_v1(s.paras)
        v2_claims = find_uncited_claims_v2(s.paras)
        v1_total += len(v1_claims)
        v2_total += len(v2_claims)
        if v2_claims:
            sub_to_claims.append({
                "subsection": s.heading,
                "claims": v2_claims,
            })

    intermediate = OUT / "_ch1_surviving_claims.json"
    with intermediate.open("w", encoding="utf-8") as fh:
        json.dump(sub_to_claims, fh, indent=2, ensure_ascii=False)
    print(f"Task 4: v1 would surface {v1_total} claims; v2 (filtered) keeps "
          f"{v2_total} -> {intermediate.name}")

    # ---- Task 5 ----
    residues = scan_numbered_residues(thesis, chapters)
    residues_md = OUT / "numbered_citation_residues.md"
    write_numbered_residues(residues_md, residues)
    chap_set = {r["chapter"] for r in residues}
    print(f"Task 5: {len(residues)} residual numbered citation(s) across "
          f"{len(chap_set)} chapter section(s) -> {residues_md.name}")

    # ---- diff summary ----
    print()
    print("=" * 72)
    print(
        f"Diff summary: the new exclusion rules removed "
        f"{v1_total - v2_total} false-positive claim(s) "
        f"(v1: {v1_total} → v2: {v2_total}). Task 5 surfaced "
        f"{len(residues)} residual non-Harvard citation(s) across "
        f"{len(chap_set)} chapter section(s). Hand-crafted topic "
        f"summaries are still required — see "
        f"{intermediate.name}."
    )


if __name__ == "__main__":
    main()
