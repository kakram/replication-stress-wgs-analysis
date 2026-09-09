"""Tests for scripts/build_annotated_references.py."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from build_annotated_references import (  # noqa: E402
    clean_crossref_abstract,
    fuzzy_threshold_for,
    gate_candidate,
    normalise_title,
    requires_surname_match,
    surnames_match,
    title_word_count,
)


class TestCleanCrossrefAbstract(unittest.TestCase):
    def test_single_entity(self):
        self.assertEqual(
            clean_crossref_abstract("foo &lt; bar"),
            "foo < bar",
        )

    def test_double_encoded_entity(self):
        """The canonical motivating case: '&amp;lt;' is the &-escaped form of
        '&lt;', which decodes to '<'. We need both passes to land at '<'."""
        self.assertEqual(
            clean_crossref_abstract("c.1988+5G&amp;lt;A"),
            "c.1988+5G<A",
        )

    def test_strips_jats_tags(self):
        self.assertEqual(
            clean_crossref_abstract(
                "<jats:p>Hello <jats:italic>world</jats:italic>.</jats:p>"
            ),
            "Hello world .",
        )

    def test_none_input(self):
        self.assertIsNone(clean_crossref_abstract(None))

    def test_empty_input(self):
        self.assertIsNone(clean_crossref_abstract(""))


class TestNormaliseTitle(unittest.TestCase):
    def test_lowercases_and_strips_punctuation(self):
        self.assertEqual(
            normalise_title("CDKN1A/p21WAF1, RB1: insights!"),
            "cdkn1a p21waf1 rb1 insights",
        )

    def test_collapses_whitespace(self):
        self.assertEqual(
            normalise_title("foo     bar\n\tbaz"),
            "foo bar baz",
        )

    def test_handles_none(self):
        self.assertEqual(normalise_title(None), "")


class TestShortTitleGating(unittest.TestCase):
    """Hardening against short, generic titles pulling in the wrong paper
    (the canonical motivating case: 'CRISPR-Cas systems' fuzzy-matching to
    a Bifidobacterium paper by Lugli instead of the intended Barrangou
    chapter)."""

    def test_word_count(self):
        self.assertEqual(title_word_count("CRISPR-Cas systems"), 2)
        self.assertEqual(title_word_count("A very long and detailed title"), 6)
        self.assertEqual(title_word_count(None), 0)
        self.assertEqual(title_word_count(""), 0)

    def test_threshold_escalates_for_short_titles(self):
        self.assertEqual(fuzzy_threshold_for("CRISPR-Cas systems"), 95)
        self.assertEqual(fuzzy_threshold_for("R loops"), 95)
        self.assertEqual(
            fuzzy_threshold_for("A very long and detailed paper title"),
            85,
        )

    def test_requires_surname_match_for_short_titles(self):
        self.assertTrue(requires_surname_match("CRISPR-Cas systems"))
        self.assertFalse(requires_surname_match(
            "A specific and reasonably long title that should be distinctive"
        ))

    def test_surnames_match_accent_insensitive(self):
        self.assertTrue(surnames_match("García", "garcia"))
        self.assertTrue(surnames_match("Álvarez-Fernández",
                                       "Alvarez-Fernandez"))
        self.assertTrue(surnames_match("Barrangou", "BARRANGOU"))
        self.assertFalse(surnames_match("Barrangou", "Lugli"))
        self.assertFalse(surnames_match("Smith", ""))
        self.assertFalse(surnames_match("", "Smith"))

    def test_short_title_rejects_surname_mismatch(self):
        """The canonical case: short title 'CRISPR-Cas systems',
        candidate authored by Lugli (Bifidobacterium paper), reference
        authored by Barrangou. Reject the match even at score 95."""
        self.assertFalse(gate_candidate(
            ref_title="CRISPR-Cas systems",
            ref_surname="Barrangou",
            cand_title="CRISPR-Cas systems in Bifidobacterium",
            cand_surname="Lugli",
            score=95,
        ))

    def test_short_title_accepts_surname_match(self):
        self.assertTrue(gate_candidate(
            ref_title="CRISPR-Cas systems",
            ref_surname="Barrangou",
            cand_title="CRISPR-Cas systems",
            cand_surname="Barrangou",
            score=95,
        ))

    def test_short_title_rejects_below_95(self):
        """Surnames match but fuzzy score is below the elevated threshold."""
        self.assertFalse(gate_candidate(
            ref_title="CRISPR-Cas systems",
            ref_surname="Barrangou",
            cand_title="CRISPR-Cas systems",
            cand_surname="Barrangou",
            score=90,
        ))

    def test_long_title_passes_at_85(self):
        """Long titles are inherently distinctive — accept at score 85
        regardless of surname (we do not enforce surname-match for them)."""
        self.assertTrue(gate_candidate(
            ref_title="A specific long descriptive paper title about "
                      "replication stress and ZFP36L1",
            ref_surname="Smith",
            cand_title="A specific long descriptive paper title about "
                       "replication stress and ZFP36L1",
            cand_surname="Smith",
            score=85,
        ))


if __name__ == "__main__":
    unittest.main()
