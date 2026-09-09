# Bibliography Merge Report

- **Timestamp:** 20260511_001751
- **Mode:** DRY-RUN (no files written)
- **references.docx:** `/Users/akramk/Documents/dev/replication-stress-wgs-analysis/thesis/references.docx`
- **referencesDelta.docx:** `/Users/akramk/Documents/dev/replication-stress-wgs-analysis/thesis/referencesDelta.docx`
- **Entries before:** 486
- **Entries after (projected):** 489
- **Shared keys (skipped):** 8
- **New entries inserted:** 3

## New entries — insertion plan

| # | Surname | Year | Original delta format | Normalised form | Insertion point | Italics preserved? | Deviations flagged |
|---|---------|------|-----------------------|-----------------|-----------------|--------------------|-------------------|
| 1 | Banerji | 2012 | `bare_year` — Banerji, S., Cibulskis, K., Rangel-Escareno, C., Brown, K.K., Carter, S.L., Frederick, A.M., Lawrence, M.S.... | Banerji, S., Cibulskis, K., Rangel-Escareno, C., Brown, K.K., Carter, S.L., Frederick, A.M., Lawrence, M.S.... | before para 25 (Barlow 2013) | no — plain text (italics lost) | — |
| 2 | Diffley | 2011 | `other` — Diffley, John FX. "Quality control in the initiation of eukaryotic DNA replication." Philosophical Transact... | Diffley, John FX. "Quality control in the initiation of eukaryotic DNA replication." Philosophical Transact... | before para 98 (Dinges 2019) | yes (runs copied) | uses double quotes (Westminster prefers single quotes); uses 'no.' for issue number (MLA-style remnant) |
| 3 | Oberle | 1991 | `bare_year` — Oberle, I., Rousseau, F., Heitz, D., Kretz, C., Devys, D., Hanauer, A., Boue, J., Bertheas, M. and Mandel, ... | Oberle, I., Rousseau, F., Heitz, D., Kretz, C., Devys, D., Hanauer, A., Boue, J., Bertheas, M. and Mandel, ... | before para 313 (Ong 2011) | no — plain text (italics lost) | — |

## Entries flagged for manual style review

| Surname | Year | Reason |
|---------|------|--------|
| Banerji | 2012 | italics lost during year-format normalisation (re-apply if needed) |
| Diffley | 2011 | year-format not auto-normalisable (needs manual rewrite to Harvard); uses double quotes (Westminster prefers single quotes); uses 'no.' for issue number (MLA-style remnant) |
| Oberle | 1991 | italics lost during year-format normalisation (re-apply if needed) |

## Entries that could not be parsed

_None._

## references.docx paragraphs the parser skipped

_5 non-empty paragraph(s) in `references.docx` could not be parsed as a reference. These are NOT modified by this script._

## Run log

```
========================================================================
merge_references.py — 2026-05-11T00:17:51
  references : /Users/akramk/Documents/dev/replication-stress-wgs-analysis/thesis/references.docx
  delta      : /Users/akramk/Documents/dev/replication-stress-wgs-analysis/thesis/referencesDelta.docx
  dry-run    : True
========================================================================
references.docx     : 492 paragraphs (492 non-empty)
referencesDelta.docx: 14 paragraphs (11 non-empty)
[dry-run] Backup would be written to references_backup_20260511_001751.docx
Parsed references   : 486 entries (unparseable paragraphs: 5)
Parsed delta        : 11 entries (unparseable paragraphs: 0)
Shared keys         : 8
NEW keys to insert  : 3
[dry-run] No file changes written.
After merge (proj.) : 489 entries
```
