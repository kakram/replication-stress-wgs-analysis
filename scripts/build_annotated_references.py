"""
Build an annotated bibliography by enriching references.docx entries with
abstracts (CrossRef → PubMed) and human-written thesis-relevance notes.

Two-phase workflow:
  Phase 1 (default):
    - Parse references.docx (uses scripts/merge_references.py parser).
    - For each entry (subject to --limit / --offset), fetch abstract from
      CrossRef; fall back to PubMed; cache results in sqlite.
    - Write outputs/annotation_coverage_report.md (the coverage table).
    - Write outputs/_pilot_annotations_intermediate.json — the per-entry data
      Claude needs in order to author relevance notes.
    - Write thesis/referencesAnnotated.docx using placeholder relevance notes.
  Phase 2 (re-run with --notes-file):
    - Read the same intermediate JSON.
    - Look up each entry's hand-written relevance note from the notes JSON.
    - Re-write thesis/referencesAnnotated.docx with real notes substituted.

Honest by design: if the abstract can't be located, the entry is emitted as
"Annotation pending — abstract not available. Reason: …" rather than
fabricated.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

from docx import Document
from rapidfuzz import fuzz

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from merge_references import (  # noqa: E402
    norm,
    parse_references_doc,
    strip_accents,
)

DEFAULT_REFS = REPO / "thesis" / "references.docx"
DEFAULT_OUT = REPO / "thesis" / "referencesAnnotated.docx"
CACHE_DIR = REPO / "scripts" / "cache"
CACHE_DB = CACHE_DIR / "abstracts.sqlite"
OUT = REPO / "outputs"
COVERAGE_REPORT = OUT / "annotation_coverage_report.md"
INTERMEDIATE_JSON = OUT / "_pilot_annotations_intermediate.json"

FUZZY_THRESHOLD = 85
SHORT_TITLE_WORD_COUNT = 4    # < this many words → "short title" handling
SHORT_TITLE_FUZZY_THRESHOLD = 95
RATE_LIMIT_SLEEP = 0.25       # CrossRef / PubMed: ~4 req/sec
SS_RATE_LIMIT_SLEEP = 0.6     # Semantic Scholar: ~1.6 req/sec (free tier)
SS_BACKOFF_429 = 5.0
RETRY_BACKOFF = 2.0
HTTP_TIMEOUT = 20


# ============================ data classes =============================== #

@dataclass
class Entry:
    paragraph_index: int
    surname: str
    year: str
    key: str
    full_text: str
    title_fragment: Optional[str] = None
    journal_fragment: Optional[str] = None
    cache_key: str = ""
    doi: Optional[str] = None
    abstract: Optional[str] = None
    source: Optional[str] = None  # "crossref" | "pubmed" | "semanticscholar" | None
    unresolved_reason: Optional[str] = None
    raw_paragraph: object = None  # Paragraph object for run-formatting copy
    cache_hit: bool = False        # set True if served from cache this run


@dataclass
class CacheStats:
    hits: int = 0
    network_fetches: int = 0
    first_time_caches: int = 0


# ============================ cache ====================================== #

def open_cache() -> sqlite3.Connection:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(CACHE_DB)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS abstracts (
            cache_key TEXT PRIMARY KEY,
            doi TEXT,
            title TEXT,
            abstract TEXT,
            source TEXT,
            fetched_at TEXT,
            raw_response TEXT
        )
    """)
    conn.commit()
    return conn


def cache_get(conn, cache_key: str) -> Optional[dict]:
    row = conn.execute(
        "SELECT doi, title, abstract, source, fetched_at, raw_response "
        "FROM abstracts WHERE cache_key = ?",
        (cache_key,),
    ).fetchone()
    if row is None:
        return None
    return {"doi": row[0], "title": row[1], "abstract": row[2],
            "source": row[3], "fetched_at": row[4],
            "raw_response": row[5]}


def cache_put(conn, cache_key: str, doi: Optional[str], title: Optional[str],
              abstract: Optional[str], source: str,
              raw_response: Optional[str]) -> None:
    conn.execute("""
        INSERT OR REPLACE INTO abstracts
            (cache_key, doi, title, abstract, source, fetched_at, raw_response)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    """, (cache_key, doi, title, abstract, source,
          datetime.now().isoformat(timespec="seconds"),
          raw_response[:8000] if raw_response else None))
    conn.commit()


# ============================ field extraction =========================== #

YEAR_PARENS_RE = re.compile(r"\(((?:19|20)\d{2}[a-z]?)\)\s*\.?\s*")
YEAR_BARE_RE = re.compile(r",\s*((?:19|20)\d{2}[a-z]?)\s*\.\s*")


def extract_title_journal(full_text: str) -> tuple[Optional[str], Optional[str]]:
    """Heuristic split: locate the year, take the next sentence as title,
    and the following segment up to the first comma as journal name."""
    t = norm(full_text)
    m = YEAR_PARENS_RE.search(t) or YEAR_BARE_RE.search(t)
    if not m:
        return None, None
    after = t[m.end():].strip()
    # Title ends at the first ". " followed by an uppercase letter — that
    # almost always begins the journal/source name.
    parts = re.split(r"\.\s+(?=[A-ZÀ-ɏ])", after, maxsplit=1)
    title = parts[0].strip(" .,'\"")
    rest = parts[1] if len(parts) > 1 else ""
    journal = rest.split(",")[0].strip() if rest else None
    if title.endswith("?") or title.endswith("!"):
        title = title  # keep punctuation
    return title or None, journal or None


def compute_cache_key(surname: str, year: str, title: Optional[str]) -> str:
    title_norm = (title or "").lower()
    title_norm = re.sub(r"[^a-z0-9]+", "", title_norm)[:80]
    blob = f"{surname.lower()}|{year}|{title_norm}"
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


# ============================ HTTP helpers =============================== #

def _ua() -> str:
    mailto = os.environ.get("CROSSREF_MAILTO", "").strip()
    base = "AnnotatedRefs/1.0 (PhD thesis tooling; +https://example.invalid)"
    if mailto:
        return f"{base} mailto:{mailto}"
    return base


def _http_get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": _ua()})
    with urllib.request.urlopen(req, timeout=HTTP_TIMEOUT) as r:
        return r.read()


def _http_get_json(url: str) -> dict:
    return json.loads(_http_get(url).decode("utf-8"))


def _retryable_call(fn, *args, **kwargs):
    """Try fn(*args, **kwargs) once, retry once after RETRY_BACKOFF on failure."""
    try:
        return fn(*args, **kwargs)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        time.sleep(RETRY_BACKOFF)
        try:
            return fn(*args, **kwargs)
        except Exception as e2:  # noqa: BLE001
            raise RuntimeError(f"network failure (after retry): {e2}") from e2


# ============================ CrossRef =================================== #

def title_word_count(title: Optional[str]) -> int:
    return len(title.split()) if title else 0


def fuzzy_threshold_for(ref_title: Optional[str]) -> int:
    """Short reference titles are inherently ambiguous; require a tighter
    fuzzy score to accept a candidate match."""
    if title_word_count(ref_title) < SHORT_TITLE_WORD_COUNT:
        return SHORT_TITLE_FUZZY_THRESHOLD
    return FUZZY_THRESHOLD


def requires_surname_match(ref_title: Optional[str]) -> bool:
    """Short titles also require an exact first-author surname match against
    the candidate, otherwise a generic title like 'CRISPR-Cas systems' can
    pull in the wrong paper."""
    return title_word_count(ref_title) < SHORT_TITLE_WORD_COUNT


def _norm_surname(s: Optional[str]) -> str:
    if not s:
        return ""
    return strip_accents(s).lower().strip()


def surnames_match(a: Optional[str], b: Optional[str]) -> bool:
    """Case- and accent-insensitive exact surname match."""
    na = _norm_surname(a)
    nb = _norm_surname(b)
    return bool(na) and na == nb


def gate_candidate(*, ref_title: Optional[str], ref_surname: str,
                   cand_title: Optional[str],
                   cand_surname: Optional[str],
                   score: int) -> bool:
    """Decide whether a fuzzy-matched candidate should be accepted.

    Long titles (>= SHORT_TITLE_WORD_COUNT words): pass with score >=
    FUZZY_THRESHOLD.
    Short titles (< 4 words): require score >= SHORT_TITLE_FUZZY_THRESHOLD
    AND that the candidate's first-author surname matches the reference's.
    """
    if score < fuzzy_threshold_for(ref_title):
        return False
    if requires_surname_match(ref_title):
        if not surnames_match(ref_surname, cand_surname):
            return False
    return True


def clean_crossref_abstract(raw: Optional[str]) -> Optional[str]:
    if not raw:
        return None
    text = re.sub(r"<[^>]+>", " ", raw)
    # CrossRef JATS payloads can carry double-encoded entities
    # (e.g. "&amp;lt;" should resolve to "<"). Unescape until stable.
    prev = None
    while text != prev:
        prev = text
        text = html.unescape(text)
    text = re.sub(r"\s+", " ", text).strip()
    return text or None


def normalise_title(s: Optional[str]) -> str:
    """Lowercase, strip punctuation, collapse whitespace — used as input to
    fuzzy-match scoring so superficial formatting differences don't matter."""
    if not s:
        return ""
    s = s.lower()
    s = re.sub(r"[^\w\s]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def crossref_lookup(entry: Entry) -> Optional[dict]:
    """Return {doi, title, abstract, raw} for the top hit if fuzzy title
    matches >= FUZZY_THRESHOLD; otherwise None."""
    q = urllib.parse.urlencode({
        "query.bibliographic": entry.full_text,
        "rows": 3,
    })
    url = f"https://api.crossref.org/works?{q}"
    raw = _retryable_call(_http_get_json, url)
    items = raw.get("message", {}).get("items", [])
    if not items:
        return None
    if not entry.title_fragment:
        # accept top hit without fuzzy gate when we have no title to match
        top = items[0]
    else:
        our_norm = normalise_title(entry.title_fragment)
        scored = []
        for item in items:
            t_list = item.get("title", [])
            title = (t_list[0] if t_list else "") or ""
            score = fuzz.token_set_ratio(our_norm, normalise_title(title))
            scored.append((score, item))
        scored.sort(key=lambda x: x[0], reverse=True)
        top_score, top = scored[0]
        cand_title_list = top.get("title", []) or []
        cand_title = cand_title_list[0] if cand_title_list else ""
        cand_authors = top.get("author", []) or []
        cand_surname = (cand_authors[0].get("family") if cand_authors else None)
        if not gate_candidate(
            ref_title=entry.title_fragment, ref_surname=entry.surname,
            cand_title=cand_title, cand_surname=cand_surname,
            score=top_score,
        ):
            return None
    doi = top.get("DOI")
    title_list = top.get("title", [])
    title = title_list[0] if title_list else None
    abstract = clean_crossref_abstract(top.get("abstract"))
    return {"doi": doi, "title": title, "abstract": abstract,
            "raw": json.dumps(top)[:6000]}


# ============================ PubMed ===================================== #

def pubmed_lookup(entry: Entry) -> Optional[dict]:
    if not entry.title_fragment:
        return None
    # PubMed's [Title] field parser silently fails on long queries
    # (empirically: ~12+ words → zero hits even for an exact match). Use the
    # first 8 distinctive words — long enough to disambiguate, short enough
    # for the parser. Surname disambiguates further.
    title_short = " ".join(entry.title_fragment.split()[:8])
    term = f'{title_short}[Title] AND {entry.surname}[Author]'
    es_q = urllib.parse.urlencode({
        "db": "pubmed", "term": term, "retmode": "json", "retmax": 3,
    })
    es_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{es_q}"
    es = _retryable_call(_http_get_json, es_url)
    idlist = es.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return None
    time.sleep(RATE_LIMIT_SLEEP)
    ef_q = urllib.parse.urlencode({
        "db": "pubmed", "id": idlist[0],
        "rettype": "abstract", "retmode": "xml",
    })
    ef_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{ef_q}"
    xml_bytes = _retryable_call(_http_get, ef_url)
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        return None
    title_el = root.find(".//ArticleTitle")
    article_title = "".join(title_el.itertext()).strip() if title_el is not None else ""
    abstract_parts = []
    for ab in root.findall(".//AbstractText"):
        label = ab.get("Label")
        text = "".join(ab.itertext()).strip()
        abstract_parts.append(f"{label}: {text}" if label else text)
    abstract = " ".join(abstract_parts).strip() or None
    score = fuzz.token_set_ratio(
        normalise_title(entry.title_fragment),
        normalise_title(article_title),
    )
    au_el = root.find(".//AuthorList/Author/LastName")
    cand_surname = au_el.text if au_el is not None else None
    if not gate_candidate(
        ref_title=entry.title_fragment, ref_surname=entry.surname,
        cand_title=article_title, cand_surname=cand_surname,
        score=score,
    ):
        return None
    pmid = idlist[0]
    return {"doi": None, "title": article_title, "abstract": abstract,
            "raw": f"PMID:{pmid}; xml_len={len(xml_bytes)}",
            "pmid": pmid}


# ============================ Semantic Scholar =========================== #

def _ss_get_json(url: str) -> dict:
    """GET that handles 429 with a single 5s backoff retry."""
    try:
        return _http_get_json(url)
    except urllib.error.HTTPError as e:
        if e.code == 429:
            time.sleep(SS_BACKOFF_429)
            return _http_get_json(url)
        raise


def semantic_scholar_lookup(entry: Entry) -> Optional[dict]:
    if not entry.title_fragment:
        return None
    q = urllib.parse.urlencode({
        "query": entry.title_fragment[:200],
        "fields": "abstract,title,authors,year,externalIds",
        "limit": 3,
    })
    url = f"https://api.semanticscholar.org/graph/v1/paper/search?{q}"
    raw = _retryable_call(_ss_get_json, url)
    data = raw.get("data", []) or []
    if not data:
        return None
    our_norm = normalise_title(entry.title_fragment)
    scored = []
    for item in data:
        title = item.get("title") or ""
        score = fuzz.token_set_ratio(our_norm, normalise_title(title))
        scored.append((score, item))
    scored.sort(key=lambda x: x[0], reverse=True)
    top_score, top = scored[0]
    cand_authors = top.get("authors", []) or []
    cand_surname = None
    if cand_authors:
        name = cand_authors[0].get("name", "") or ""
        parts = name.split()
        # Semantic Scholar author names are typically "FirstName LastName" —
        # take the last whitespace-separated token as a surname proxy.
        cand_surname = parts[-1] if parts else None
    if not gate_candidate(
        ref_title=entry.title_fragment, ref_surname=entry.surname,
        cand_title=top.get("title"), cand_surname=cand_surname,
        score=top_score,
    ):
        return None
    abstract = (top.get("abstract") or "").strip() or None
    if abstract is None:
        return None  # title matches but no abstract — no point caching
    doi = (top.get("externalIds") or {}).get("DOI")
    return {"doi": doi, "title": top.get("title"), "abstract": abstract,
            "raw": json.dumps(top)[:6000]}


# ============================ fetch orchestration ======================== #

def _recover_cached_reason(raw_response: Optional[str]) -> Optional[str]:
    """Pull a previously-stored failure reason out of raw_response JSON, if any."""
    if not raw_response:
        return None
    try:
        data = json.loads(raw_response)
    except (json.JSONDecodeError, TypeError):
        return None
    if isinstance(data, dict) and "reason" in data:
        return data["reason"]
    return None


def fetch_one(entry: Entry, conn, *,
              cache_only: bool, no_cache: bool,
              stats: CacheStats, ss_enabled: bool) -> None:
    """Populate entry.{abstract, source, doi, unresolved_reason}."""
    cached = None if no_cache else cache_get(conn, entry.cache_key)
    if cached and cached.get("abstract"):
        entry.abstract = cached["abstract"]
        entry.doi = cached.get("doi")
        entry.source = cached.get("source") or "cache"
        entry.cache_hit = True
        stats.hits += 1
        return
    if cache_only:
        # Prefer recovering the original first-pass reason from raw_response.
        if cached:
            recovered = _recover_cached_reason(cached.get("raw_response"))
            entry.unresolved_reason = recovered or (
                "no cached abstract; --cache-only set"
            )
            entry.doi = cached.get("doi")
        else:
            entry.unresolved_reason = "no cached abstract; --cache-only set"
        return

    stats.network_fetches += 1
    if not cached:
        stats.first_time_caches += 1

    # ---- CrossRef ----
    try:
        cr = crossref_lookup(entry)
        cr_err = None
    except Exception as e:  # noqa: BLE001
        cr = None
        cr_err = str(e)
    time.sleep(RATE_LIMIT_SLEEP)

    if cr and cr.get("abstract"):
        entry.abstract = cr["abstract"]
        entry.doi = cr.get("doi")
        entry.source = "crossref"
        cache_put(conn, entry.cache_key, cr.get("doi"), cr.get("title"),
                  cr.get("abstract"), "crossref", cr.get("raw"))
        return

    # ---- PubMed ----
    try:
        pm = pubmed_lookup(entry)
        pm_err = None
    except Exception as e:  # noqa: BLE001
        pm = None
        pm_err = str(e)
    time.sleep(RATE_LIMIT_SLEEP)

    if pm and pm.get("abstract"):
        entry.abstract = pm["abstract"]
        entry.doi = cr.get("doi") if cr else None
        entry.source = "pubmed"
        cache_put(conn, entry.cache_key,
                  entry.doi, pm.get("title"), pm.get("abstract"),
                  "pubmed", pm.get("raw"))
        return

    # ---- Semantic Scholar (3rd tier, gated by SS_ENABLED) ----
    if ss_enabled:
        try:
            ss = semantic_scholar_lookup(entry)
            ss_err = None
        except Exception as e:  # noqa: BLE001
            ss = None
            ss_err = str(e)
        time.sleep(SS_RATE_LIMIT_SLEEP)
    else:
        ss = None
        ss_err = None

    if ss and ss.get("abstract"):
        entry.abstract = ss["abstract"]
        entry.doi = ss.get("doi") or (cr.get("doi") if cr else None)
        entry.source = "semanticscholar"
        cache_put(conn, entry.cache_key, entry.doi, ss.get("title"),
                  ss.get("abstract"), "semanticscholar", ss.get("raw"))
        return

    # ---- exhausted; record the most informative failure reason ----
    reasons = []
    if cr is None:
        reasons.append(f"crossref error: {cr_err[:80]}" if cr_err
                       else "crossref: no fuzzy-match hit")
    elif not cr.get("abstract"):
        reasons.append("crossref: matched DOI but no abstract field")
    if pm is None:
        reasons.append(f"pubmed error: {pm_err[:80]}" if pm_err
                       else "pubmed: no fuzzy-match hit")
    elif not pm.get("abstract"):
        reasons.append("pubmed: matched but abstract empty")
    if not ss_enabled:
        reasons.append("semanticscholar: disabled (SS_ENABLED env var unset)")
    elif ss is None:
        reasons.append(f"semanticscholar error: {ss_err[:80]}" if ss_err
                       else "semanticscholar: no fuzzy-match hit")
    elif not ss.get("abstract"):
        reasons.append("semanticscholar: matched but abstract empty")
    reason = "; ".join(reasons) or "no hits"
    entry.unresolved_reason = reason

    # If CrossRef recorded a DOI, hold onto it for the docx output.
    entry.doi = cr.get("doi") if cr else None

    # Cache the failure with the reason so --cache-only can recover it later.
    cache_put(conn, entry.cache_key, entry.doi, cr.get("title") if cr else None,
              None, "unresolved", json.dumps({"reason": reason}))


# ============================ docx output ================================ #

PLACEHOLDER_NOTE = "[RELEVANCE NOTE TO BE WRITTEN]"


def copy_runs(source_p, target_p):
    for run in source_p.runs:
        nr = target_p.add_run(run.text)
        if run.italic is not None:
            nr.italic = run.italic
        if run.bold is not None:
            nr.bold = run.bold
        if run.underline is not None:
            nr.underline = run.underline


SOURCE_TAGS = {
    "crossref":          "[CrossRef]",
    "pubmed":            "[PubMed]",
    "semanticscholar":   "[Semantic Scholar]",
    "cache":             "[Cached]",
    "manual_pdf_extract": "[Manual (PDF)]",
}


def _note_text(note_data) -> str:
    """Notes JSON values may be a legacy plain string (the placeholder) or
    the current dict shape {note, source, abstract, ...}. Extract the note
    text robustly from either."""
    if isinstance(note_data, dict):
        return note_data.get("note") or PLACEHOLDER_NOTE
    return note_data or PLACEHOLDER_NOTE


def _note_abstract(note_data) -> Optional[str]:
    if isinstance(note_data, dict):
        return note_data.get("abstract")
    return None


def _note_source(note_data) -> Optional[str]:
    if isinstance(note_data, dict):
        return note_data.get("source")
    return None


def render_annotation_text(entry: Entry, note_data) -> str:
    note_text = _note_text(note_data)
    note_abstract = _note_abstract(note_data)
    note_source = _note_source(note_data)

    # If the entry is unresolved but the notes file supplies a manually-
    # extracted abstract (e.g. predecessor theses), prefer that abstract +
    # note over the generic "Annotation pending" stub.
    if entry.unresolved_reason and not note_abstract:
        return ("Annotation pending — abstract not available. "
                f"Reason: {entry.unresolved_reason}")

    if entry.abstract:
        abstract = entry.abstract
        source_tag = SOURCE_TAGS.get(entry.source or "", "[?]")
        doi_str = f" doi:{entry.doi}" if entry.doi else ""
    else:
        abstract = note_abstract or "(empty abstract)"
        source_tag = SOURCE_TAGS.get(note_source or "", "[Manual]")
        doi_str = ""

    return (f"{source_tag}{doi_str} Abstract: {abstract}\n\n"
            f"Relevance: {note_text}")


def build_annotated_docx(entries: list[Entry], notes: dict[str, str],
                         output_path: Path) -> None:
    doc = Document()
    doc.add_heading("Annotated References", level=0)
    doc.add_paragraph(
        "Each entry is followed by the abstract retrieved from CrossRef or "
        "PubMed (attributed in square brackets), and a thesis-relevance "
        "note. Entries for which an abstract could not be located are "
        "marked 'Annotation pending'."
    )
    for entry in entries:
        # Reference paragraph (preserve runs/italic where possible).
        ref_p = doc.add_paragraph()
        if entry.raw_paragraph is not None:
            copy_runs(entry.raw_paragraph, ref_p)
        else:
            ref_p.add_run(entry.full_text)
        # Annotation paragraph.
        note = notes.get(entry.cache_key, PLACEHOLDER_NOTE)
        ann_text = render_annotation_text(entry, note)
        ann_p = doc.add_paragraph()
        ann_p.add_run(ann_text).italic = False
        # Spacer.
        doc.add_paragraph("")
    doc.save(str(output_path))


# ============================ reports / intermediate ===================== #

def _normalised_reason(reason: Optional[str]) -> str:
    """Reduce a verbose semicolon-separated reason to a coarse bucket
    so the halt-summary can group similar failures together."""
    if not reason:
        return "unknown"
    if "crossref: matched DOI but no abstract" in reason and \
            "pubmed: no fuzzy-match" in reason:
        return "CrossRef DOI matched but no abstract field; PubMed no hit"
    if "no fuzzy-match" in reason and "no fuzzy-match" in reason:
        return "CrossRef and PubMed: no fuzzy-match hit"
    if "matched but abstract empty" in reason:
        return "matched but abstract empty"
    if "error" in reason:
        return "network/API error"
    return reason.split(";")[0].strip()


def _print_halt_summary(entries: list[Entry]) -> None:
    by_src: dict[str, int] = {}
    reasons: dict[str, int] = {}
    for e in entries:
        if e.abstract:
            by_src[e.source or "unknown"] = by_src.get(e.source or "unknown", 0) + 1
        else:
            by_src["unresolved"] = by_src.get("unresolved", 0) + 1
            r = _normalised_reason(e.unresolved_reason)
            reasons[r] = reasons.get(r, 0) + 1
    print()
    print("Per-source breakdown:")
    for src in ("crossref", "pubmed", "semanticscholar", "unresolved"):
        print(f"  {src:<18} {by_src.get(src, 0)}")
    print()
    print("Top 5 reasons for unresolved:")
    top = sorted(reasons.items(), key=lambda x: x[1], reverse=True)[:5]
    for reason, n in top:
        print(f"  {n:>4}  {reason}")


def write_coverage_report(entries: list[Entry], stats: CacheStats,
                          path: Path) -> None:
    by_src = {"crossref": 0, "pubmed": 0, "semanticscholar": 0, "unresolved": 0}
    for e in entries:
        if e.abstract:
            key = e.source if e.source in by_src else "unresolved"
            by_src[key] += 1
        else:
            by_src["unresolved"] += 1
    total = len(entries) or 1
    seconds_saved = stats.hits * 0.5  # rough placeholder

    with path.open("w", encoding="utf-8") as fh:
        fh.write("# Annotation Coverage Report\n\n")
        fh.write(f"- **Generated:** "
                 f"{datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"- **Entries processed:** {len(entries)}\n\n")

        fh.write("## Source breakdown (where did the data originate)\n\n")
        fh.write("| Source           | Count | %    |\n")
        fh.write("|------------------|------:|-----:|\n")
        for key, label in [
            ("crossref",         "CrossRef"),
            ("pubmed",           "PubMed"),
            ("semanticscholar",  "Semantic Scholar"),
            ("unresolved",       "Unresolved"),
        ]:
            n = by_src.get(key, 0)
            pct = 100 * n / total
            fh.write(f"| {label:<16} | {n:>5} | {pct:>4.1f} |\n")
        fh.write("\n")

        fh.write("## Cache utilisation (this run)\n\n")
        fh.write(f"- Cache hits: **{stats.hits}** "
                 f"(saved approx. {seconds_saved:.0f} seconds at "
                 f"0.5 s/entry network cost)\n")
        fh.write(f"- Network fetches: **{stats.network_fetches}**\n")
        fh.write(f"- First-time entries: **{stats.first_time_caches}** "
                 "(cached for next run)\n\n")

        unresolved = [e for e in entries if not e.abstract]
        fh.write("## Unresolved entries\n\n")
        if not unresolved:
            fh.write("_None._\n")
        else:
            fh.write("| Surname | Year | Reason |\n")
            fh.write("|---------|------|--------|\n")
            for e in unresolved:
                reason = (e.unresolved_reason or "").replace("|", "\\|")
                fh.write(f"| {e.surname} | {e.year} | {reason} |\n")


def write_intermediate_json(entries: list[Entry], path: Path) -> None:
    payload = []
    for e in entries:
        payload.append({
            "cache_key": e.cache_key,
            "surname": e.surname,
            "year": e.year,
            "key": e.key,
            "title_fragment": e.title_fragment,
            "journal_fragment": e.journal_fragment,
            "doi": e.doi,
            "source": e.source,
            "abstract": e.abstract,
            "unresolved_reason": e.unresolved_reason,
            "full_text_first_180": e.full_text[:180],
        })
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False),
                    encoding="utf-8")


# ============================ main pipeline ============================== #

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--references", type=Path, default=DEFAULT_REFS)
    p.add_argument("--output", type=Path, default=DEFAULT_OUT)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--offset", type=int, default=0)
    p.add_argument("--cache-only", action="store_true")
    p.add_argument("--no-cache", action="store_true")
    p.add_argument("--notes-file", type=Path, default=None,
                   help="JSON file mapping cache_key → relevance note text.")
    p.add_argument("--allow-low-coverage", action="store_true",
                   help="Skip the 70%% coverage halt-check; write the docx "
                        "even if coverage is below threshold. Use only after "
                        "reviewing the halt report.")
    args = p.parse_args(argv)

    if args.cache_only and args.no_cache:
        print("FATAL: --cache-only and --no-cache are mutually exclusive.")
        return 2

    if not args.references.exists():
        print(f"FATAL: references file not found at {args.references}")
        return 2

    ss_enabled = os.environ.get("SS_ENABLED", "").strip() in (
        "1", "true", "yes", "True", "YES"
    )
    if ss_enabled:
        print("Semantic Scholar fallback: ENABLED")
    else:
        print("Semantic Scholar fallback: DISABLED "
              "(set SS_ENABLED=1 to enable)")

    refs_doc = Document(str(args.references))
    parsed, unparseable = parse_references_doc(refs_doc)
    print(f"Parsed {len(parsed)} entries from {args.references.name} "
          f"({len(unparseable)} unparseable, skipped).")

    # Apply --offset and --limit.
    slice_ = parsed[args.offset:]
    if args.limit is not None:
        slice_ = slice_[: args.limit]
    print(f"Processing slice: offset={args.offset}, limit={args.limit}, "
          f"count={len(slice_)}.")

    # Build Entry objects with title/journal extracted.
    entries: list[Entry] = []
    for ref in slice_:
        full_text = ref.text
        title, journal = extract_title_journal(full_text)
        entry = Entry(
            paragraph_index=ref.paragraph_index,
            surname=ref.surname,
            year=ref.year,
            key=ref.key,
            full_text=full_text,
            title_fragment=title,
            journal_fragment=journal,
            cache_key=compute_cache_key(ref.surname, ref.year, title),
            raw_paragraph=ref.paragraph,
        )
        entries.append(entry)

    # Fetch.
    conn = open_cache()
    stats = CacheStats()
    for i, e in enumerate(entries, 1):
        fetch_one(e, conn,
                  cache_only=args.cache_only, no_cache=args.no_cache,
                  stats=stats, ss_enabled=ss_enabled)
        status = e.source or "unresolved"
        tag = " (cache)" if e.cache_hit else ""
        print(f"  [{i:3d}/{len(entries)}] {e.surname:<22} {e.year}  "
              f"-> {status}{tag}"
              + (f"  ({len(e.abstract)} chars)" if e.abstract else ""))

    # ---- 70% coverage halt-check ----
    resolved = sum(1 for e in entries if e.abstract)
    total = len(entries) or 1
    coverage_pct = 100.0 * resolved / total
    print()
    print(f"Coverage: {resolved}/{total} ({coverage_pct:.1f}%).")

    # Always write the coverage report + intermediate JSON so the user can
    # diagnose; only the docx is gated.
    write_coverage_report(entries, stats, COVERAGE_REPORT)
    write_intermediate_json(entries, INTERMEDIATE_JSON)

    if coverage_pct < 70.0 and not args.allow_low_coverage:
        print()
        print("=" * 72)
        print(f"HALT: coverage {coverage_pct:.1f}% is below the 70% threshold.")
        print("referencesAnnotated.docx was NOT written.")
        print(f"See {COVERAGE_REPORT.relative_to(REPO)} for full diagnostics.")
        _print_halt_summary(entries)
        print()
        print("To override (e.g. after manual review), re-run with "
              "--allow-low-coverage.")
        print("=" * 72)
        return 4

    # Load relevance notes if provided.
    notes: dict[str, str] = {}
    if args.notes_file and args.notes_file.exists():
        notes = json.loads(args.notes_file.read_text(encoding="utf-8"))
        print(f"Loaded {len(notes)} relevance note(s) from "
              f"{args.notes_file.name}.")

    build_annotated_docx(entries, notes, args.output)

    print()
    print(f"Coverage report: {COVERAGE_REPORT.relative_to(REPO)}")
    print(f"Intermediate JSON: {INTERMEDIATE_JSON.relative_to(REPO)}")
    print(f"Annotated docx: {args.output.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
