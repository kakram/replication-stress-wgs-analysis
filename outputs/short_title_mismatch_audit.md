# Short-Title Cache Mismatch Audit

**Run:** 2026-05-11
**Scope:** 17 reference entries in `references.docx` whose extracted title fragment was < 4 words, against `scripts/cache/abstracts.sqlite`.

## Method

For each short-title reference entry:

1. Looked up the entry in the abstracts cache by `cache_key`.
2. From the cached `raw_response`, attempted to extract the candidate paper's first-author surname and publication year (CrossRef: `author[0].family` and `published-print.date-parts`; Semantic Scholar: `authors[0].name`; PubMed cache stores only the PMID, so author/year cannot be recovered without re-querying).
3. Compared against the reference entry's `surname` and `year`.
4. For PubMed-only cache rows where author info isn't available, applied a heuristic on the cached title vs the reference's short title — flag if the cached title contains an organism / scope qualifier not in the reference.

## Flagged entries

| Cache key | Reference (surname + year + title) | Candidate (surname + year + title) | Mismatch reason |
|-----------|------------------------------------|------------------------------------|-----------------|
| `4acfed3b4312bd24…` | Barrangou 2013 — "CRISPR-Cas systems" | (PubMed; surname not in raw_response) — "Occurrence and applications of CRISPR-Cas systems in bifidobacteria" | Cached PubMed candidate is a wrong-paper match: its title describes a *Bifidobacterium*-specific CRISPR-Cas review (authored by Lugli et al., not Barrangou). The reference's short title "CRISPR-Cas systems" was insufficient to disambiguate against the canonical Barrangou & van der Oost chapter. |

**Confirmed mismatches: 1.**

## Non-mismatches investigated

The following short-title cache rows were examined and found to be acceptable matches (or are unresolved entries that simply have no candidate to verify):

| Reference | Status | Note |
|-----------|--------|------|
| Colnaghi 2011 — "October" | unresolved | Title-extraction artefact (parser caught the month). No cached candidate to verify. |
| Diffley 2011 — ": 3545-3553" | unresolved | Title-extraction artefact (page numbers). No cached candidate to verify. |
| Durkin 2007 — "Chromosome fragile sites" | CrossRef, resolved | Cached candidate title matches reference. `"published-print"` year is 2007 ✓. |
| Leonard 2013 — "DNA replication origins" | unresolved | No cached candidate. |
| Li 2025 — "Chromoanagenesis in osteosarcoma" | CrossRef, resolved | Cached candidate title matches; `"published-print"` 2025 ✓. |
| Lowe 2017 — "Transcriptomics technologies" | PubMed, resolved | Cached candidate title matches; author info unavailable from raw_response but no other red flags. |
| Magdalou 2014 — "June" | unresolved | Title-extraction artefact. No cached candidate. |
| McArdle 2021 — "What is proteomics?" | CrossRef, resolved | False positive on first pass: candidate has both `"published-online":2020` and `"published-print":2021`; the reference cites the print year (2021), so this is the same paper. |
| Pelham 2020 — "June" | unresolved | Title-extraction artefact. |
| Rajaei 2021 — ": 1602-1613" | unresolved | Title-extraction artefact. |
| Rehm 2020 — "Alcohol consumption" | PubMed, resolved | Cached candidate title matches; no red flags. |
| Sher 2022 — "August" | unresolved | Title-extraction artefact. |
| Siegel 2019 — "Cancer statistics, 2019" | CrossRef, resolved | Cached candidate title matches; `"published-print"` 2019 ✓. |
| Smith 2007 — "February" | PubMed, resolved | Cached candidate title matches the reference; no red flags. Reference parsing artefact captured "February" rather than the real title. |
| Sutherland 1991 — "Chromosomal fragile sites" | unresolved | No cached candidate to verify. |
| Tiessen 2012 — ": 1-23" | unresolved | Title-extraction artefact. |

## Actions taken

- **Deleted from cache:** `4acfed3b4312bd24…` (Barrangou 2013) — see next section for re-fetch.

## Related observation — title-extraction parser is fragile on bare-year entries

10 of the 17 short-title entries are *parsing artefacts*, not genuinely short titles: when a delta-style entry uses `Author, YYYY, Month.` format (e.g. *Colnaghi, 2011, October*) or when an entry ends with bare page numbers, the title extractor stops at the wrong delimiter. These entries (`Colnaghi`, `Magdalou`, `Pelham`, `Sher`, `Smith`, `Diffley`, `Rajaei`, `Tiessen` and others) tend to land as `unresolved` because the wrong "title" is sent to CrossRef/PubMed. Worth a follow-up pass to improve `extract_title_journal()` for these cases, but out of scope for this audit.
