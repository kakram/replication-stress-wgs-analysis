# §5.10 SIFT del:tol Fisher exact tests — Summary

Source: `outputs/ch5_variant_burden.csv` (sift_deleterious, sift_tolerated).

## SIFT counts per sample

| Sample | SIFT-del | SIFT-tol | Total | % del |
|---|---:|---:|---:|---:|
| WTUN | 213 | 381 | 594 | 35.9% |
| WTAPH | 197 | 385 | 582 | 33.8% |
| B4UN | 36 | 35 | 71 | 50.7% |
| B4APH | 32 | 29 | 61 | 52.5% |
| **WT pooled** | **410** | **766** | **1176** | **34.9%** |
| **B4 pooled** | **68** | **64** | **132** | **51.5%** |

Bonferroni-corrected threshold for the three pairwise tests: αₜₒᵣᵣ = 0.05 / 3 = **0.0167**.

## Test results

### Test 1: WTUN vs WTAPH (within-WT treatment)

- Group 1: 213 / 594 (35.9%) deleterious
- Group 2: 197 / 582 (33.8%) deleterious
- Odds ratio: **1.092**  (95% CI: 0.853, 1.399)
- p-value: **0.5009** — not significant

### Test 2: B4UN vs B4APH (within-B4 treatment)

- Group 1: 36 / 71 (50.7%) deleterious
- Group 2: 32 / 61 (52.5%) deleterious
- Odds ratio: **0.933**  (95% CI: 0.444, 1.955)
- p-value: **0.8629** — not significant

### Test 3: WT pooled vs B4 pooled (between-genotype)

- Group 1: 410 / 1176 (34.9%) deleterious
- Group 2: 68 / 132 (51.5%) deleterious
- Odds ratio: **0.504**  (95% CI: 0.345, 0.736)
- p-value: **0.0002655** — **SIGNIFICANT**

### Test 4: 2x4 chi-square across all four samples

- χ² = 14.745, dof = 3, p = 0.002048
