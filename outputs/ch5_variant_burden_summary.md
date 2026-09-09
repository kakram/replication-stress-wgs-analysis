# §5.3 MCF-7 Somatic Variant Burden — Summary

Source data: Sentieon TNhaplotyper2 + TNfilter somatic calls for MCF-7 WT (parental) and B4 (CRISPR/Cas9 ZFP36L1$^{-/-}$ clone, from Teotia 2024), each $\pm$ aphidicolin. Files used per sample: `<SAMPLE>_somatic.vcf.gz` (FILTER classes), `<SAMPLE>_somatic_vep_anno.maf.gz` (coding / SIFT / gnomAD), `<SAMPLE>_tmb.tsv` (TMB).

## Numbers

| Metric | WTUN | WTAPH | B4UN | B4APH |
|---|---|---|---|---|
| Total variants (VCF) | 1,198,693 | 1,153,726 | 249,482 | 212,078 |
| PASS variants | 144,145 | 141,954 | 19,463 | 18,326 |
|   FILTER: germline | 509,461 | 489,428 | 26,801 | 27,322 |
|   FILTER: panel_of_normals | 0 | 0 | 0 | 0 |
|   FILTER: weak_evidence | 173,415 | 174,906 | 96,761 | 85,782 |
|   FILTER: clustered_events | 574,099 | 549,416 | 116,894 | 99,164 |
|   FILTER: other | 583,633 | 554,734 | 211,385 | 169,234 |
| Coding variants (total) | 1,143 | 1,139 | 144 | 130 |
|   Missense | 613 | 600 | 71 | 61 |
|   Synonymous | 331 | 334 | 44 | 40 |
|   Nonsense | 35 | 35 | 7 | 9 |
|   Frameshift deletion | 12 | 10 | 0 | 1 |
|   Frameshift insertion | 10 | 11 | 1 | 0 |
|   In-frame deletion | 7 | 5 | 3 | 0 |
|   In-frame insertion | 9 | 8 | 0 | 2 |
|   Splice site | 13 | 10 | 4 | 10 |
|   Splice region | 111 | 124 | 14 | 7 |
| SIFT deleterious | 213 | 197 | 36 | 32 |
| SIFT tolerated | 381 | 385 | 35 | 29 |
| SIFT unknown | 138,217 | 135,992 | 19,060 | 17,889 |
| gnomAD known (AF > 0.001) | 315 | 323 | 25 | 37 |
| gnomAD novel (AF ≤ 0.001) | 138,496 | 136,251 | 19,106 | 17,913 |
| TMB (mut/Mb) | 0.42 | 0.35 | 0.07 | 0.05 |
| Genes with ≥1 variant | 21,080 | 20,972 | 6,738 | 6,334 |
| Genes with coding variant | 815 | 816 | 129 | 116 |

## WT-vs-B4 ratios

```
WT-U vs B4-U:
  PASS variants ratio   :  144,145 /   19,463  = 7.41
  Coding variants ratio :    1,143 /      144  = 7.94
  TMB ratio             :     0.42 /     0.07  = 6.00

WT-A vs B4-A:
  PASS variants ratio   :  141,954 /   18,326  = 7.75
  Coding variants ratio :    1,139 /      130  = 8.76
  TMB ratio             :     0.35 /     0.05  = 7.00
```

## Per-sample summary

```
Sample          PASS    Coding   Missense   SIFT-del     TMB
WTUN         144,145     1,143        613        213    0.42
WTAPH        141,954     1,139        600        197    0.35
B4UN          19,463       144         71         36    0.07
B4APH         18,326       130         61         32    0.05
```

## Diagnostic

The WT-vs-B4 asymmetry is **consistent across treatment**: the PASS-variant ratio is 7.41 untreated vs 7.75 aphidicolin-treated (within ±20 %), and the TMB ratio shows the same pattern (6.00 vs 7.00). This is the signature of a **baseline-burden / clonal-heterogeneity effect** rather than a genotype × treatment interaction: the WT and B4 lines differ in their somatic-variant load irrespective of aphidicolin exposure. Plausible causes include the WT line being a polyclonal MCF-7 stock against which the B4 single-cell-derived ZFP36L1$^{-/-}$ clone is a bottlenecked subline (so 'WT > B4' may simply reflect pre-existing clonal diversity), or systematic differences in coverage / variant-calling sensitivity between the two line states. The §5.3 prose should foreground this clonality caveat before claiming any aphidicolin or ZFP36L1-loss effect from these counts.
