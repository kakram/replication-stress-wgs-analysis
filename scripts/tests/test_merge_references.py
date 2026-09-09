"""Minimal tests for scripts/merge_references.py.

Run with: python -m unittest scripts.tests.test_merge_references
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

from docx import Document

# Make the script importable when running directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from merge_references import (  # noqa: E402
    detect_year_format,
    normalise_to_parens_year,
    parse_ref_paragraph,
    run,
)


def _make_refs_docx(path: Path, entries: list[str]) -> None:
    doc = Document()
    doc.add_paragraph("References")
    for e in entries:
        doc.add_paragraph(e)
    doc.save(str(path))


class TestYearFormat(unittest.TestCase):
    def test_parens_form_is_canonical(self):
        text = "Adams, J. (2018). A paper. Journal, 1(1), pp.1-10."
        self.assertEqual(detect_year_format(text), "parens")

    def test_bare_year_with_month(self):
        text = "Chase, T., 2020, January. A paper. Journal, 1(1), pp.1-10."
        self.assertEqual(detect_year_format(text), "bare_year_month")

    def test_bare_year_with_period(self):
        text = "Davis, M., 2019. A paper. Journal, 2(2), pp.20-30."
        self.assertEqual(detect_year_format(text), "bare_year")

    def test_mla_style_is_other(self):
        text = ('Diffley, John FX. "Quality control in DNA replication." '
                'Phil. Trans. R. Soc. B 366, no. 1584 (2011): 3545-3553.')
        self.assertEqual(detect_year_format(text), "other")


class TestNormalisation(unittest.TestCase):
    def test_bare_year_month_normalisation(self):
        text = "Chase, T., 2020, January. Topic. Journal, 1(1), pp.1-10."
        out = normalise_to_parens_year(text, "2020", "bare_year_month")
        self.assertIn("Chase, T. (2020).", out)
        self.assertNotIn("January", out)

    def test_bare_year_normalisation(self):
        text = "Davis, M., 2019. Topic. Journal, 2(2), pp.20-30."
        out = normalise_to_parens_year(text, "2019", "bare_year")
        self.assertIn("Davis, M. (2019).", out)

    def test_parens_unchanged(self):
        text = "Adams, J. (2018). Topic. Journal, 1(1), pp.1-10."
        out = normalise_to_parens_year(text, "2018", "parens")
        self.assertEqual(out, text)


class TestParseRefParagraph(unittest.TestCase):
    def test_parens_form(self):
        self.assertEqual(
            parse_ref_paragraph("Adams, J. (2018). Topic. Journal."),
            ("Adams", "2018"),
        )

    def test_bare_year_with_month(self):
        self.assertEqual(
            parse_ref_paragraph("Chase, T., 2020, January. Topic. Journal."),
            ("Chase", "2020"),
        )

    def test_diacritic_led_surname(self):
        # Á (U+00C1) at start of surname must match.
        self.assertEqual(
            parse_ref_paragraph(
                "Álvarez-Fernández, M., Sanz-Flores, M. and Manchado, E. "
                "(2017). Therapeutic relevance. Cancer Cell, 31(5), pp.669-684."
            ),
            ("Álvarez-Fernández", "2017"),
        )

    def test_diacritic_in_middle_of_surname(self):
        # Lüönd (U+0308 combining diaeresis) — diacritics inside surname,
        # ASCII first letter — already covered by \w; this guards against
        # regressions in the inside-token character class.
        self.assertEqual(
            parse_ref_paragraph(
                "Lüönd, F., Tiede, S. and Christofori, G. (2021). "
                "Breast cancer as an example. British Journal of Cancer, 124, "
                "pp.164-175."
            ),
            ("Lüönd", "2021"),
        )

    def test_year_volume_not_confused_with_year(self):
        # Nature volume 550, issue 7676 — must not be mistaken for year.
        self.assertEqual(
            parse_ref_paragraph(
                "Shendure, J. (2017). DNA sequencing at 40. "
                "Nature, 550(7676), pp.345-353."
            ),
            ("Shendure", "2017"),
        )


class TestDryRunInsertion(unittest.TestCase):
    """End-to-end dry-run: 5 ref entries + 2 delta entries; verify insertion
    points and year-format normalisation."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.tmp.name)
        self.refs_path = self.tmp_path / "references.docx"
        self.delta_path = self.tmp_path / "referencesDelta.docx"
        self.report_dir = self.tmp_path / "outputs"
        self.report_dir.mkdir()

        _make_refs_docx(self.refs_path, [
            "Adams, J. (2018). Paper one. Journal, 1(1), pp.1-10.",
            "Brown, K. (2019). Paper two. Journal, 2(2), pp.20-30.",
            "Carter, L. (2017). Paper three. Journal, 3(3), pp.40-50.",
            "Davis, M. (2020). Paper four. Journal, 4(4), pp.60-70.",
            "Evans, N. (2021). Paper five. Journal, 5(5), pp.80-90.",
        ])
        _make_refs_docx(self.delta_path, [
            # Already conformant; should slot between Adams and Brown.
            "Black, S. (2020). Paper six. Journal, 6(6), pp.100-110.",
            # Bare-year-month; should normalise + slot between Carter and Davis.
            "Chase, T., 2020, January. Paper seven. Journal, 7(7), "
            "pp.120-130.",
        ])

    def tearDown(self):
        self.tmp.cleanup()

    def test_dry_run_inserts_alphabetically_and_normalises(self):
        rc = run(
            references_path=self.refs_path,
            delta_path=self.delta_path,
            dry_run=True,
            report_dir=self.report_dir,
        )
        self.assertEqual(rc, 0, "merge run should exit zero in dry-run")

        # references.docx must be untouched in dry-run.
        doc_after = Document(str(self.refs_path))
        body = [p.text for p in doc_after.paragraphs if p.text.strip()]
        self.assertEqual(len(body), 6, "dry-run must not modify the file")
        self.assertEqual(body[0], "References")

        # Read the report and verify content.
        reports = list(self.report_dir.glob("merge_references_report_*.md"))
        self.assertEqual(len(reports), 1, "exactly one report should be written")
        text = reports[0].read_text(encoding="utf-8")

        # Both new entries listed.
        self.assertIn("Black", text)
        self.assertIn("Chase", text)

        # Year-format detection.
        self.assertIn("`parens`", text)
        self.assertIn("`bare_year_month`", text)

        # Year-format normalisation applied for Chase.
        self.assertIn("Chase, T. (2020).", text)

        # Alphabetical insertion targets: Black before Brown, Chase before Davis.
        # The "References" header occupies paragraph 0, so:
        #   Adams=1, Brown=2, Carter=3, Davis=4, Evans=5.
        self.assertIn("before para 2 (Brown 2019)", text)
        self.assertIn("before para 4 (Davis 2020)", text)


if __name__ == "__main__":
    unittest.main()
