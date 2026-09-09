# PhD Thesis Project — ZFP36L1 WGS Analysis

## Author & context
PhD thesis by Khalid Akram (University of Westminster). Investigating the
genomic consequences of ZFP36L1 loss in MCF-7 (breast cancer) and U2OS
(osteosarcoma) cell lines under aphidicolin-induced replication stress.
WGS performed by Azenta/GeneWiz.

Predecessor theses in the lab (cite in thesis):
- Solaiman (2021) — U2OS ZFP36L1-/- generation, cytology
- Sidali (2023) — R-loops, FRA14C, chromatin association
- Teotia (2024) — MCF-7 B4 clone generation

## Current focus: Chapter 6 — U2OS WGS analysis
MCF-7 analysis (Chapter 5) is complete and lives in this repo.
U2OS analysis (Chapter 6) is starting now. Strategy: mirror the MCF-7
analysis pipeline so Chapter 7 (cross-lineage comparison) becomes natural.

## Samples
- MCF-7: WT-UN, WT-APH, B4-UN, B4-APH (B4 = ZFP36L1-/- clone, from Teotia)
- U2OS:  WT-U,  WT-A,  G3-U,  G3-A   (G3 = ZFP36L1-/- clone, from Solaiman)

## Data layout — IMPORTANT
The full U2OS Azenta deliverable lives on an external drive and is NOT
inside this repo. Only the files needed for active analysis are copied
into data/u2os/ on demand. Do not assume any specific U2OS file is
present — check first, and if missing, ask before doing anything else.

Repo data layout:
  data/u2os/                  <- subset of U2OS files copied from external drive
    <files copied as needed, e.g.:>
    WT-U_vep_anno_non_common.tsv.gz
    WT-A_vep_anno_non_common.tsv.gz
    G3-U_vep_anno_non_common.tsv.gz
    G3-A_vep_anno_non_common.tsv.gz
    <and more added as analyses require>

  data/                       <- existing MCF-7 data (pre-existing)
  data/cfs_extracted/         <- common fragile site BED files
  data/CFS_fixed/

  results/u2os/               <- U2OS analysis outputs (create as needed)
  results/                    <- existing MCF-7 outputs

  figures/                    <- final figures for thesis
  scripts/                    <- analysis scripts (existing MCF-7 + new U2OS)

When a script needs a U2OS file that isn't in data/u2os/ yet, STOP and
tell me what filename(s) to bring over from the external drive. Do not
assume paths outside the repo.

## Methodological note (IMPORTANT — flag in Methods chapter)
MCF-7 was processed with Sentieon TNseq (somatic caller).
U2OS was processed with germline caller + non_common (rare-variant) filter.
Pipelines are not identical. Document in Methods, acknowledge in Discussion.
Frame Chapter 7 cross-lineage comparison as "trends" not exact-match numbers.

## Established U2OS findings so far (from VEP all_stats.json)
Knockout (G3) vs wild-type effect on the non_common filtered set:
- Total filtered records:   +1.6%   (WT mean 803,390  -> G3 mean 816,466)
- Missense variants:        +16.9%  (1,850 -> 2,163)
- Frameshift variants:      +50.7%  (104 -> 156)   [HEADLINE FINDING]
- HIGH impact variants:     +8.5%   (441 -> 478)
- MODERATE impact:          +14.0%  (2,438 -> 2,778)
- SIFT-deleterious:         +22.4%  (584 -> 715)
Aphidicolin effect within each genotype is small (a few hundred variants).
Genotype is the dominant axis of variation, not drug treatment.

Chromosomal redistribution observed at unfiltered level (likely karyotypic):
chr3 +10.1%, chr15 +13.9%, chr18 -15.6%. Needs CNV check to confirm.

## Immediate todo (in order)
1. Top-30 mutated genes per U2OS sample (from non_common.tsv.gz)
   Required files: data/u2os/{WT-U,WT-A,G3-U,G3-A}_vep_anno_non_common.tsv.gz
2. Combined side-by-side top-10 table for Chapter 6 Section 6.4
3. Rainfall plots (4 U2OS conditions) — reuse/parametrise variantvis.R
   Required files: same as above (positions extractable from same TSVs)
4. CNV analysis — small files, easy win
   Required files: data/u2os/{WT-U,WT-A,G3-U,G3-A}_cnv.cns
                   data/u2os/{WT-U,WT-A,G3-U,G3-A}_cnv.vcf.gz
5. CFS overlap for U2OS top genes (reuse compute_cfs_overlaps.py)
6. SV analysis from data/u2os/*_sv.vcf.gz

## Existing reusable scripts in this repo (verify before reusing)
- scripts/compute_cfs_overlaps.py    — CFS overlap analysis (MCF-7)
- scripts/compute_genes_near_CFS.py  — gene-CFS proximity
- scripts/heatmap_mutated_genes.py   — heatmap visualisation
- scripts/upset_plot_mutated_genes.py — UpSet plot
- scripts/plot_cfs_ideogram_hg38_scaled.py
- scripts/plot_genomewide_ideogram.py
- variantvis.R                       — R rainfall plot script

When working on U2OS, parametrise these existing scripts where possible
rather than duplicating. Goal: same code, different inputs.

## Conventions
- Reference genome: GRCh38/hg38
- U2OS outputs go to results/u2os/, MCF-7 outputs stay in results/
- Figures go to figures/ (subfolders u2os/ and mcf7/ if helpful)
- Python: pandas, numpy, matplotlib, pysam if needed
- R: tidyverse, Bioconductor (karyoploteR, GenomicRanges)
- Commit working analyses; don't batch unrelated changes
- Show diffs before applying changes to existing files

## Style for code-related questions
- Always check what U2OS files actually exist in data/u2os/ before
  writing code that depends on them. If a needed file is missing,
  STOP and tell me what to bring over.
- Prefer reusing/parametrising existing MCF-7 scripts over writing new ones
- When writing new code, match the style of existing scripts in this repo
- Always print sanity-check counts (rows in, rows out, genes found, etc.)

## Bibliography Management — Invariants

These rules apply to ALL operations that modify references.docx. They are non-negotiable.

1. **Format.** Every entry in references.docx MUST follow Westminster-Harvard style 
   (parenthetical-year Harvard, per the University of Westminster Harvard Referencing 
   Quick Guide, copy at: Harvard_Referencing_Quick_Guide_2025_9_1.pdf in project root). 
   The canonical year format is `(YYYY)`, NOT `, YYYY.` or `, YYYY, Month.`

2. **Order.** Entries MUST be inserted in alphabetical order by first-author surname, 
   then by year ascending within the same surname. NEVER append new entries at the 
   end of the file.

3. **Canonical tool.** ALL merges from referencesDelta.docx into references.docx MUST 
   go through `scripts/merge_references.py`. Do not write inline merge logic in 
   one-off scripts or sessions. If the canonical script doesn't handle a case, 
   extend the script rather than working around it.

4. **Safety.** The canonical script must back up references.docx (SHA-256 verified) 
   before any modification and must refuse to proceed if backup verification fails.

5. **Semantic Scholar fallback** is OFF by default due to aggressive rate-limiting on the free tier. Re-enable with SS_ENABLED=1 only when an S2 API key is available (set S2_API_KEY env var).

6. **Relevance notes** can be generated in two modes: `--mode interactive` (CC session writes them) or `--mode api` (Sonnet 4.6 via Anthropic SDK). Use api mode for bulk generation; interactive for hand-crafted notes on high-stakes references like predecessor theses.

## Deferred — IGV deep-dives for Chapter 6 (after MCF-7 backfill complete)
Genes to interrogate in IGV across all four U2OS BAMs (WT-U, WT-A, G3-U, G3-A):
- FHIT (chr3:59.7-61.2 Mb, FRA3B)         — headline finding
- ZFP36L1 (chr14:68.79 Mb, FRA14C)        — knockout validation
- FAM8A1 (chr6:17.6 Mb)                   — strongest coding candidate
- FCGBP (chr19:39.9 Mb, FRA19A)           — single-replicate mystery
- ROBO2 (chr3:75.9-77.7 Mb)               — chr3 corroborating finding

Output: 200-300 words + 2-4 IGV screenshots per gene = §6.6 Gene-level deep dives
Format: parallel to Chapter 5's CLCN3/HRNR/MT-ND5/ZNF180/FAM186A treatment
- Use relative paths from repo root (e.g. data/u2os/..., not absolute paths)
