"""
Relevance-note generator for the annotated bibliography.

This is the persistence/orchestration layer. The actual *writing* of each
60-100-word relevance note is done by Claude Code (the assistant), reading
the abstract for each entry and producing a thesis-themed note. The script:

  1. Loads outputs/_pilot_annotations_intermediate.json (abstracts).
  2. Loads outputs/relevance_notes.json (existing notes).
  3. Computes the set of entries needing a note for this batch
     (--offset / --limit; resolved entries only unless overridden).
  4. Writes outputs/_relevance_notes_todo_batch.json — the work queue
     Claude consumes.
  5. After Claude writes outputs/relevance_notes.json, re-running with
     --verify reports gaps still remaining in the batch.

Schema for outputs/relevance_notes.json:

  {
    "<cache_key>": {
        "surname": "...", "year": "...",
        "note": "60-100 word relevance note",
        "source": "claude_from_abstract" | "manual_pdf_extract" | ...,
        "abstract": "...",                 # optional override; used for
                                           # predecessor theses where the
                                           # cache has no abstract
        "notes_generated_at": "2026-05-11T..."
    },
    ...
  }
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

REPO = Path(__file__).resolve().parent.parent
DEFAULT_INTERMEDIATE = REPO / "outputs" / "_pilot_annotations_intermediate.json"
DEFAULT_NOTES = REPO / "outputs" / "relevance_notes.json"
DEFAULT_TODO_BATCH = REPO / "outputs" / "_relevance_notes_todo_batch.json"

FLUSH_EVERY = 25

# ---- API mode constants ----
API_MODEL = "claude-sonnet-4-6"
API_MAX_TOKENS = 300
API_TEMPERATURE = 0.4
API_RATE_LIMIT_SLEEP = 1.5  # <50 req/min → ≥1.2s; 1.5s is comfortable headroom
API_RETRY_BACKOFF = 5.0
# Sonnet 4.6 pricing (per million tokens)
PRICE_INPUT = 3.0
PRICE_OUTPUT = 15.0
PRICE_CACHE_WRITE = 3.75   # 1.25× base
PRICE_CACHE_READ = 0.30    # 0.10× base

# Few-shot exemplar cache keys (stable across runs). Drawn from the existing
# notes.json — one foundational, one partial-relevance, one background.
FEW_SHOT_KEYS = {
    "foundational":      "ac14138a1a66c23dab20825c77f071e4a66aa92aa1b8a9ccfc44ee6538b11fd4",  # Bartkova 2005
    "partial-relevance": "7334c5b0a253f0211e1a8c9f76666c0c4fdf75db08140d59697c1f81a37797ed",  # Aftab 2021
    "background":        "1d641a29ecfca2a6de6f04ac5b405d74b54b144e656c43dcafd9844fa7d1b0c1",  # Burguin 2021
}


# ---------------------------- IO helpers ---------------------------------- #

def load_intermediate(path: Path) -> list[dict]:
    if not path.exists():
        raise SystemExit(f"FATAL: intermediate JSON not found at {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_notes(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    raw = json.loads(path.read_text(encoding="utf-8"))
    # Tolerate legacy str-keyed-by-cache_key {cache_key: "note_text"} format.
    if raw and isinstance(next(iter(raw.values())), str):
        raw = {k: {"note": v, "source": "legacy", "notes_generated_at": ""}
               for k, v in raw.items()}
    return raw


def save_notes(notes: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(notes, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    tmp.replace(path)


# ---------------------------- selection ----------------------------------- #

def select_batch(intermediate: list[dict], notes: dict[str, dict],
                 offset: int, limit: int | None,
                 regenerate: bool) -> list[dict]:
    """Return entries needing a note for this batch, in alphabetical order.

    Default behaviour: only entries with abstracts (resolved) are included.
    Predecessor theses (no abstract) are handled separately by --include-unresolved.
    """
    resolved = [e for e in intermediate if e.get("abstract")]
    resolved.sort(key=lambda e: (e["surname"].lower(), e["year"]))
    sliced = resolved[offset:]
    if limit is not None:
        sliced = sliced[:limit]
    if regenerate:
        return sliced
    return [e for e in sliced if e["cache_key"] not in notes]


# ============================ API mode =================================== #

def confirm_api_key() -> str:
    """Confirm ANTHROPIC_API_KEY is set in env or .env at repo root."""
    # Try loading .env if present, without forcing it.
    env_path = REPO / ".env"
    if env_path.exists():
        try:
            from dotenv import load_dotenv  # type: ignore
            load_dotenv(env_path)
        except ImportError:
            pass
    key = os.environ.get("ANTHROPIC_API_KEY", "").strip()
    if not key:
        raise SystemExit(
            "FATAL: ANTHROPIC_API_KEY not found in environment or "
            f"{env_path}. Halting before any API calls. "
            "Set the key (export ANTHROPIC_API_KEY=sk-...) and re-run."
        )
    print(f"ANTHROPIC_API_KEY: set (len={len(key)}).")
    return key


THESIS_CONTEXT_PREAMBLE = """\
You write thesis-relevance notes for a PhD bibliography. The PhD investigates
the genomic consequences of ZFP36L1 loss in MCF-7 (breast cancer) and U2OS
(osteosarcoma) cell lines under aphidicolin-induced replication stress, using
whole-genome sequencing (WGS).

Each note connects the supplied paper's abstract to one or more of these
thesis themes, naming the relevant chapter / subsection where useful:

  Themes:
    - ZFP36L1 / RNA-binding protein biology
    - Replication stress, aphidicolin, DNA damage response (DDR)
    - Common fragile sites (CFSs), FRA3B, FHIT, FRA14C
    - MCF-7 / U2OS / breast cancer / osteosarcoma biology
    - Whole-genome sequencing, mutational burden, variant interpretation
    - CRISPR/Cas9 genome editing methodology
    - Clonal evolution, tumour heterogeneity

  Chapter map:
    Ch1 Introduction — has H2 subsections numbered 1.1-1.15 covering
        hallmarks of cancer (1.1), breast cancer (1.2), animal / 2D / 3D
        models (1.3-1.4), CRISPR history & applications (1.5-1.7),
        BC + DDR (1.8), tumour evolution (1.9), ZFP36L1 (1.10), cell
        plasticity (1.11), DNA replication & replication stress (1.12),
        aphidicolin (1.13), WGS + genome instability (1.14), WGS in
        cancer research (1.15).
    Ch5 — MCF-7 WGS deep-dives, including candidate-gene subsections for
        CLCN3, HRNR, MT-ND5, ZNF180, FAM186A.
    Ch6 — U2OS WGS, with CNV and CFS-overlap subsections; G3 clone is
        the ZFP36L1-/- system from Solaiman (2021).
    Ch7 — Comparative MCF-7 vs U2OS synthesis.
    Ch8 — General discussion.

Honesty rules — non-negotiable:
  - If the paper is foundational to a theme, say so clearly and cite a
    specific chapter / subsection where it should be used.
  - If the paper is partially relevant (e.g. touches on breast cancer but
    not on ZFP36L1 or replication stress), say so plainly.
  - If the paper is tangential or off-theme, say so. Recommend "background
    citation only" or omission. Do not pretend tangential papers are
    central.
  - Do not invent details that aren't in the abstract.

Format:
  - 60-100 words. Single paragraph. UK English.
  - No filler ("This paper is interesting...", "This study is important...").
  - Refer to specific chapter / subsection numbers when relevant
    (e.g. "Ch1 §1.12", "Ch5 §5.4.2 CLCN3").
"""


def build_few_shot_examples(notes: dict, intermediate: list[dict]) -> str:
    """Compose the 3 input→output exemplars as a single text block."""
    inter_by_key = {e["cache_key"]: e for e in intermediate}
    blocks: list[str] = []
    for label, key in FEW_SHOT_KEYS.items():
        n = notes.get(key)
        e = inter_by_key.get(key)
        if not n or not e:
            raise SystemExit(
                f"FATAL: few-shot exemplar missing for {label} "
                f"(cache_key {key[:12]}…). Check relevance_notes.json + "
                "intermediate JSON."
            )
        blocks.append(
            f"## Example — {label}\n"
            f"Surname: {e['surname']}\n"
            f"Year: {e['year']}\n"
            f"Title: {e.get('title_fragment') or '(unknown)'}\n"
            f"Journal: {e.get('journal_fragment') or '(unknown)'}\n"
            f"Abstract: {e['abstract']}\n\n"
            f"Relevance note ({len(n['note'].split())} words):\n{n['note']}"
        )
    return "\n\n".join(blocks)


def build_user_message_content(preamble_plus_examples: str,
                               entry: dict) -> list[dict]:
    """Two-block user message: a cached prefix (preamble + examples), then
    the per-entry tail that varies and is not cached."""
    entry_block = (
        "## Current entry — write the relevance note for this one\n"
        f"Surname: {entry['surname']}\n"
        f"Year: {entry['year']}\n"
        f"Title: {entry.get('title_fragment') or '(unknown)'}\n"
        f"Journal: {entry.get('journal_fragment') or '(unknown)'}\n"
        f"Abstract: {entry['abstract']}\n\n"
        "Write the relevance note. Output the note text only, no preamble, "
        "no quotes."
    )
    return [
        {"type": "text",
         "text": preamble_plus_examples,
         "cache_control": {"type": "ephemeral"}},
        {"type": "text", "text": entry_block},
    ]


def call_api(client, entry: dict, preamble_plus_examples: str):
    """One API call with one-shot retry on transient errors."""
    import anthropic  # local import to avoid hard dep when --mode interactive

    def _do_call():
        return client.messages.create(
            model=API_MODEL,
            max_tokens=API_MAX_TOKENS,
            temperature=API_TEMPERATURE,
            system="You write factual, honest, specific relevance notes for "
                   "a PhD bibliography. Match the style and depth of the "
                   "examples provided.",
            messages=[{
                "role": "user",
                "content": build_user_message_content(
                    preamble_plus_examples, entry
                ),
            }],
        )

    try:
        resp = _do_call()
    except (anthropic.AuthenticationError,
            anthropic.PermissionDeniedError) as e:
        raise SystemExit(f"FATAL: auth/quota error — halting batch: {e}")
    except (anthropic.APIConnectionError, anthropic.APITimeoutError,
            anthropic.InternalServerError, anthropic.RateLimitError) as e:
        time.sleep(API_RETRY_BACKOFF)
        try:
            resp = _do_call()
        except (anthropic.AuthenticationError,
                anthropic.PermissionDeniedError) as ee:
            raise SystemExit(f"FATAL: auth/quota error — halting batch: {ee}")
        except anthropic.APIError as ee:
            return None, None, f"persistent API error: {ee}"

    # Successful response → extract text
    text = "".join(
        block.text for block in resp.content if getattr(block, "type", "") == "text"
    ).strip()
    return text, resp.usage, None


def generate_via_api(batch: list[dict], notes: dict, notes_path: Path,
                     intermediate: list[dict]) -> dict:
    """Drive the API-mode batch loop. Updates `notes` in place and returns
    a usage / cost summary dict."""
    import anthropic

    confirm_api_key()
    client = anthropic.Anthropic()

    preamble = THESIS_CONTEXT_PREAMBLE
    examples = build_few_shot_examples(notes, intermediate)
    preamble_plus_examples = preamble + "\n\n## Few-shot examples\n\n" + examples

    totals = {
        "input_tokens": 0,
        "cache_creation_tokens": 0,
        "cache_read_tokens": 0,
        "output_tokens": 0,
        "calls_succeeded": 0,
        "calls_failed": 0,
        "failures": [],
    }
    t0 = time.time()
    for i, entry in enumerate(batch, 1):
        text, usage, err = call_api(client, entry, preamble_plus_examples)
        if err is not None or text is None:
            totals["calls_failed"] += 1
            totals["failures"].append({"surname": entry["surname"],
                                        "year": entry["year"],
                                        "reason": err or "empty response"})
            print(f"  [{i:>3}/{len(batch)}] {entry['surname']:<22} "
                  f"{entry['year']}  -> FAIL  ({err})")
            time.sleep(API_RATE_LIMIT_SLEEP)
            continue

        wc = len(text.split())
        bucket_count = (
            usage.input_tokens + (usage.cache_creation_input_tokens or 0)
            + (usage.cache_read_input_tokens or 0)
        )
        totals["input_tokens"] += usage.input_tokens
        totals["cache_creation_tokens"] += usage.cache_creation_input_tokens or 0
        totals["cache_read_tokens"] += usage.cache_read_input_tokens or 0
        totals["output_tokens"] += usage.output_tokens
        totals["calls_succeeded"] += 1

        notes[entry["cache_key"]] = {
            "surname": entry["surname"],
            "year": entry["year"],
            "note": text,
            "source": "api_sonnet_4_6",
            "notes_generated_at": datetime.now().isoformat(timespec="seconds"),
            "tokens_in": usage.input_tokens,
            "tokens_out": usage.output_tokens,
            "cache_creation_tokens": usage.cache_creation_input_tokens or 0,
            "cache_read_tokens": usage.cache_read_input_tokens or 0,
        }
        elapsed = time.time() - t0
        print(f"  [{i:>3}/{len(batch)}] {entry['surname']:<22} "
              f"{entry['year']}  -> {wc}w  "
              f"in={usage.input_tokens}+cr={usage.cache_read_input_tokens or 0}"
              f" out={usage.output_tokens}  t={elapsed:.1f}s")

        # Incremental flush every FLUSH_EVERY successful calls.
        if totals["calls_succeeded"] % FLUSH_EVERY == 0:
            save_notes(notes, notes_path)
            print(f"  [flush] saved {len(notes)} notes to {notes_path.name}")

        time.sleep(API_RATE_LIMIT_SLEEP)

    save_notes(notes, notes_path)
    return totals


def report_cost(totals: dict) -> None:
    """Print a cost summary using two formulae: the simple one the user
    asked for ($3/M input, $15/M output, ignoring cache savings) and the
    actual billed cost using Anthropic's published cache rates."""
    raw_input_total = (
        totals["input_tokens"] + totals["cache_creation_tokens"]
        + totals["cache_read_tokens"]
    )
    simple_cost = (raw_input_total / 1e6) * PRICE_INPUT \
        + (totals["output_tokens"] / 1e6) * PRICE_OUTPUT
    actual_cost = (
        (totals["input_tokens"] / 1e6) * PRICE_INPUT
        + (totals["cache_creation_tokens"] / 1e6) * PRICE_CACHE_WRITE
        + (totals["cache_read_tokens"] / 1e6) * PRICE_CACHE_READ
        + (totals["output_tokens"] / 1e6) * PRICE_OUTPUT
    )
    savings = simple_cost - actual_cost
    print()
    print(f"API usage summary:")
    print(f"  Calls succeeded     : {totals['calls_succeeded']}")
    print(f"  Calls failed        : {totals['calls_failed']}")
    print(f"  Input tokens        : {totals['input_tokens']:>10,}")
    print(f"  Cache creation      : {totals['cache_creation_tokens']:>10,} "
          f"(written once per cache miss; reused thereafter)")
    print(f"  Cache reads         : {totals['cache_read_tokens']:>10,} "
          f"(billed at 10% of base)")
    print(f"  Output tokens       : {totals['output_tokens']:>10,}")
    print(f"  Cost (simple, $3/M input + $15/M output, ignoring cache):"
          f"     ${simple_cost:.4f}")
    print(f"  Cost (actual, incl. cache_write @ 1.25× and cache_read @ 0.1×):"
          f"    ${actual_cost:.4f}")
    if totals["cache_read_tokens"]:
        print(f"  Caching saved        : ${savings:.4f} "
              f"({100*savings/simple_cost:.1f}% of simple cost)")


def select_unresolved(intermediate: list[dict],
                      surnames: list[tuple[str, str]]) -> list[dict]:
    """Return unresolved entries matching (surname, year) tuples — used by
    --include-unresolved for predecessor thesis seeding."""
    pairs = {(s.lower(), y) for s, y in surnames}
    return [e for e in intermediate
            if not e.get("abstract")
            and (e["surname"].lower(), e["year"]) in pairs]


# ---------------------------- main ---------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--intermediate", type=Path, default=DEFAULT_INTERMEDIATE)
    p.add_argument("--notes-file", type=Path, default=DEFAULT_NOTES)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--offset", type=int, default=0)
    p.add_argument("--regenerate", action="store_true",
                   help="Re-do entries that already have notes.")
    p.add_argument("--verify", action="store_true",
                   help="Report whether the requested batch has been filled "
                        "in outputs/relevance_notes.json; do not write a "
                        "todo file.")
    p.add_argument("--todo-out", type=Path, default=DEFAULT_TODO_BATCH,
                   help="Where to write the todo-batch file consumed by "
                        "Claude.")
    p.add_argument("--mode", choices=["interactive", "api"],
                   default="interactive",
                   help="interactive (default): emit a todo file for Claude "
                        "Code to consume in-session; api: call the "
                        f"{API_MODEL} Anthropic Messages API for each entry.")
    args = p.parse_args(argv)

    intermediate = load_intermediate(args.intermediate)
    notes = load_notes(args.notes_file)
    batch = select_batch(intermediate, notes,
                         offset=args.offset, limit=args.limit,
                         regenerate=args.regenerate)

    resolved_total = sum(1 for e in intermediate if e.get("abstract"))
    print(f"Intermediate: {len(intermediate)} entries "
          f"({resolved_total} resolved, with abstracts).")
    print(f"Existing notes: {len(notes)}.")
    print(f"Batch window: offset={args.offset}, "
          f"limit={args.limit}, regenerate={args.regenerate}.")
    print(f"Entries needing a note in this batch: {len(batch)}.")

    if args.mode == "api":
        if not batch:
            print("Nothing to generate in this batch.")
            return 0
        totals = generate_via_api(batch, notes, args.notes_file, intermediate)
        report_cost(totals)
        if totals["failures"]:
            print()
            print("Failed entries:")
            for f in totals["failures"]:
                print(f"  {f['surname']} {f['year']}: {f['reason']}")
        return 0 if not totals["failures"] else 1

    if args.verify:
        # Verification mode — recompute the batch *as if* notes were empty,
        # then check how many of those entries have been filled.
        full_batch = select_batch(intermediate, {}, offset=args.offset,
                                  limit=args.limit, regenerate=False)
        filled = sum(1 for e in full_batch if e["cache_key"] in notes)
        missing = [e for e in full_batch if e["cache_key"] not in notes]
        print()
        print(f"Verify: {filled}/{len(full_batch)} notes present for this "
              "batch.")
        if missing:
            print(f"Still missing ({len(missing)}):")
            for e in missing[:20]:
                print(f"  {e['surname']} {e['year']}  ({e['cache_key'][:12]})")
            if len(missing) > 20:
                print(f"  ... and {len(missing) - 20} more")
        return 0 if not missing else 1

    # Non-verify mode: write the todo-batch JSON Claude will consume.
    todo_payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "intermediate_path": str(args.intermediate),
        "notes_path": str(args.notes_file),
        "offset": args.offset,
        "limit": args.limit,
        "regenerate": args.regenerate,
        "themes": [
            "ZFP36L1 / RNA-binding protein biology",
            "Replication stress, aphidicolin, DNA damage response",
            "Common fragile sites, FRA3B, FHIT, FRA14C",
            "MCF-7 / U2OS / breast cancer / osteosarcoma biology",
            "Whole-genome sequencing, mutational burden, variant interpretation",
            "CRISPR/Cas9 genome editing methodology",
            "Clonal evolution, tumour heterogeneity",
        ],
        "note_length_words": "60-100",
        "incremental_flush_every": FLUSH_EVERY,
        "entries": [
            {
                "cache_key": e["cache_key"],
                "surname": e["surname"],
                "year": e["year"],
                "title_fragment": e.get("title_fragment"),
                "journal_fragment": e.get("journal_fragment"),
                "doi": e.get("doi"),
                "source": e.get("source"),
                "abstract": e["abstract"],
            }
            for e in batch
        ],
    }
    args.todo_out.parent.mkdir(parents=True, exist_ok=True)
    args.todo_out.write_text(
        json.dumps(todo_payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    try:
        rel = args.todo_out.relative_to(REPO)
    except ValueError:
        rel = args.todo_out
    print(f"Todo batch written: {rel}")
    print()
    print("Next: read the abstracts in that file, write 60-100 word notes,")
    print(f"and save them to {args.notes_file.name} (dict keyed by cache_key)")
    print("with fields {note, source, abstract?, notes_generated_at}.")
    print(f"Flush every {FLUSH_EVERY} entries to avoid losing work on crash.")
    print()
    print("After saving notes, re-run with --verify to confirm completeness.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
