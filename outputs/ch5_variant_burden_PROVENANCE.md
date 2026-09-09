# Provenance — `outputs/ch5_variant_burden.csv`

Traced and verified 15 June 2026. Every figure in Chapter 5 §5.3 (and §5.3.2) derives from this CSV. This document records the full lineage from the CSV back to the Azenta source files, how each column is computed, and one methodological issue surfaced during the trace.

## 1. Lineage chain (CSV → source)

```
Azenta WGS deliverable (project 40-1022750518, hg38/hs38DH)
        │   Sentieon 202308.02: alignment → BQSR → TNhaplotyper2 → TNfilter → VEP → MAF
        ▼
data/mcf7/<SAMPLE>/
    <SAMPLE>_somatic.vcf.gz          (FILTER-class counts)
    <SAMPLE>_somatic_vep_anno.maf.gz (coding consequences, SIFT, gnomAD, gene sets)
    <SAMPLE>_tmb.tsv                  (pre-computed TMB)
        │   parsed by scripts/lib/somatic_burden.py
        │   driven by  scripts/ch5_variant_burden.py
        ▼
outputs/ch5_variant_burden.csv          ← the analysis table
outputs/ch5_variant_burden_thesis_table.tex  (Table 5.3.1)
outputs/ch5_variant_burden_summary.md
        │   read by hand
        ▼
Chapter 5 §5.3.1, §5.3.2 prose + Table 5.3.2
```

`<SAMPLE>` ∈ {WTUN, WTAPH, B4UN, B4APH}.

## 2. The four source files per sample

All four sample directories under `data/mcf7/` carry the same three files, all dated 11 May 2026 (date copied from the external drive, not the Azenta processing date).

| File | Size (WTUN) | Role | Read by |
|---|---|---|---|
| `WTUN_somatic.vcf.gz` | 65 MB | Full Mutect2/TNfilter FILTER vocabulary (germline, panel_of_normals, weak_evidence, clustered_events). The MAF collapses FILTER to PASS/common_variant, so the VCF is the only source for filter-class counts. | `parse_vcf_filters()` |
| `WTUN_somatic_vep_anno.maf.gz` | 8.2 MB | VEP-annotated MAF: `Variant_Classification`, `SIFT`, `gnomAD_AF`, `Hugo_Symbol`. | `parse_maf()` |
| `WTUN_tmb.tsv` | 53 B | Sentieon TMB: sample, target_size, count, TMB(mut/Mb). | `parse_tmb()` |

`.tbi` files are tabix indices for the VCFs (not read by the analysis).

## 3. Azenta pipeline, from the VCF headers

The `##SentieonCommandLine` headers inside each VCF record the exact processing. Verbatim, for WTUN:

- **Caller:** `TNhaplotyper2`, `sentieon-genomics-202308.02`, dated 2024-07-24.
- **Filter:** `TNfilter` (same version), using orientation-bias and contamination models.
- **Reference:** `/…/hg38/hs38DH/genome/genome.fa` (GRCh38 + decoy, hs38DH).
- **Germline resource:** `af-only-gnomad.hg38.pass.vcf.gz` (GATK gnomAD AF resource).
- **Input:** `WTUN.aln.bam` + `WTUN_BQSR.table`.
- **Project:** `40-1022750518`.

This confirms the MCF-7 pipeline stated in `CLAUDE.md` (Sentieon TNseq, hg38/hs38DH).

## 4. How each CSV column is computed

From `scripts/lib/somatic_burden.py`:

- **total_variants, pass_variants, filter_\*** — streamed from the VCF FILTER column (field 7). A multi-filter row increments *each* named bucket once (plus `filter_other` for any unnamed token), so filter columns can sum to more than `total − PASS` by design.
- **coding_total, missense, synonymous (= MAF "Silent"), nonsense, frameshift_\*, in_frame_\*, splice_\*** — counted from the MAF `Variant_Classification` column. "Coding" = the 11 classes in `CODING_CLASSES`.
- **sift_deleterious / tolerated / unknown** — MAF `SIFT` prefix match ("deleterious…" / "tolerated…" / else unknown). Non-missense rows have no SIFT score and fall into `unknown`, which is why `sift_unknown` ≈ total MAF rows.
- **gnomad_known / novel** — MAF `gnomAD_AF`; AF > 0.001 = known, else novel (`GNOMAD_NOVEL_THRESHOLD`).
- **genes_with_any_variant / coding_variant** — unique `Hugo_Symbol` sets.
- **tmb_mut_per_mb** — read directly from `<SAMPLE>_tmb.tsv` (target_size 978,354,808 bp for all four; counts 409/346/66/49).

## 5. Reproducibility — verified

Regenerated on 15 June 2026 by running `python3 scripts/ch5_variant_burden.py` against `data/mcf7/`. The output was **byte-identical** to the committed `outputs/ch5_variant_burden.csv` (`diff` reported no differences). The pipeline is pure-Python (gzip + csv), no external bioinformatics tools required to reproduce the table from the Azenta files.

```
python3 scripts/ch5_variant_burden.py
diff outputs/ch5_variant_burden.csv <committed copy>   # identical
```

## 6. ⚠ Methodological issue surfaced during the trace — the two genotypes were called differently

The `##SentieonCommandLine.TNhaplotyper2` headers and the `#CHROM` sample columns show that the four samples were **not** processed with an identical comparison scheme:

| Sample | tumor_sample | normal_sample | VCF sample cols | MAF matched-norm | Meaning of its variants |
|---|---|---|---|---|---|
| WTUN | WTUN | *(none)* | WTUN | NORMAL (placeholder) | tumour-only: variants vs hg38 reference, germline-subtracted by gnomAD AF |
| WTAPH | WTAPH | *(none)* | WTAPH | NORMAL (placeholder) | tumour-only: vs hg38 reference |
| B4UN | B4UN | **WTUN** | B4UN, WTUN | WTUN | tumour–normal: variants in B4UN **not present in WTUN** |
| B4APH | B4APH | **WTAPH** | B4APH, WTAPH | WTAPH | tumour–normal: variants in B4APH **not present in WTAPH** |

This means the WT and B4 counts are **not the same kind of measurement**. The WT call set is "everything that differs from the reference genome"; the B4 call set is "what is private to B4 after subtracting its matched WT parent". Because the B4 clone is *derived from* MCF-7 WT, a B4-vs-WT paired call removes by construction the large pool of variants the two share — which is most of them.

**Consequence for the thesis.** The 6–8× WT-vs-B4 burden gap reported in §5.3.1 is therefore driven substantially by this asymmetric calling design, not solely by the polyclonal-vs-monoclonal clonal architecture currently given as the explanation. The two explanations are entangled and the present data cannot fully separate them. This needs to be:

1. **Documented in Methods (Chapter 3)** — state explicitly that WT was called tumour-only and B4 tumour–normal against matched WT, with the Sentieon command lines as evidence.
2. **Reframed in §5.3.1** — the absolute cross-genotype burden gap is partly a pipeline artefact; cross-genotype *absolute* counts are not directly comparable. Within-genotype contrasts (WT-UN vs WT-APH; B4-UN vs B4-APH) and *proportional* analyses remain valid.
3. **Flagged in §5.3.2** — the coding-consequence *proportions* are more robust to this than the absolute counts, but the populations compared are still different in kind (B4-private somatic vs WT-vs-reference).
4. **Carried into the Chapter 7 parity discussion** alongside the existing MCF-7 (somatic) vs U2OS (germline) caveat.

This resolves the open question previously logged ("which sample served as normal in the MCF-7 TNseq pairings"): **WTUN/WTAPH had no matched normal (tumour-only); B4UN/B4APH used WTUN/WTAPH respectively as the matched normal.**

## 7. Version-control gap

As of the trace, the following are **untracked** in git (`git ls-files` returns nothing for them):
`data/mcf7/`, `scripts/ch5_variant_burden.py`, `scripts/lib/somatic_burden.py`, `outputs/ch5_variant_burden.csv`.

The analysis is reproducible only because the inputs happen to be present on disk. Committing the scripts (and recording the source-file SHA-256s, since the large `data/` files may stay out of git) would make the §5.3 chain auditable independent of the working tree.
