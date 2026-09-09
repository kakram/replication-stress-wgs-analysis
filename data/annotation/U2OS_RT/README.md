# U2OS replication-timing track — drop files here

Rename each file to the exact name below on arrival. The extension can stay
whatever it came as (.bedGraph, .bg, .bw, .txt, .gz) — the analysis script
detects the format.

## Primary: 4D Nucleome, experiment set 4DNES3Y8V593 (GRCh38)

| Save as              | What it is                          |
|----------------------|-------------------------------------|
| U2OS_early_rep1.*    | Early-S, biological replicate 1     |
| U2OS_early_rep2.*    | Early-S, biological replicate 2     |
| U2OS_late_rep1.*     | Late-S, biological replicate 1      |
| U2OS_late_rep2.*     | Late-S, biological replicate 2      |

File accessions given in the handover note (VERIFY each on the portal before
citing — the assay, cell line, build and replicate must match):

  Early-S rep1   4DNFIHZLUMUG
  Early-S rep2   4DNFIVY7ZP8F
  Late-S  rep1   4DNFIKUUEZCI
  Late-S  rep2   4DNFIXU4FWWY

Prefer the bedGraph-style processed signal over BigWig where both are offered.
Record for each file: accession, assembly, bin size, normalisation, release date.

## Optional validation: GEO GSE211592 (hg19, ratio already computed)

Three WT replicates, Early/Late ratio precomputed. hg19, but the analysis
pipeline already works in hg19 space, so these need no extra handling.

  GSE211592_RepSeq_U2OS_IDH2_W_rep1.hg19_win_EarlyLateRatio_RT.bedgraph.gz
  GSE211592_RepSeq_U2OS_IDH2_W_rep2.hg19_win_EarlyLateRatio_RT.bedgraph.gz
  GSE211592_RepSeq_U2OS_IDH2_W_rep3.hg19_win_EarlyLateRatio_RT.bedgraph.gz

Save as U2OS_hg19_ratio_rep1/2/3.bedgraph.gz
