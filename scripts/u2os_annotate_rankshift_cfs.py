#!/usr/bin/env python3
"""
Three post-processing steps on the U2OS top-30 mutated gene results.

Task 1: Annotate top30 CSVs with a 'flag' column (PAR_artefact / HLA_locus).
Task 2: Build rank_shift_u2os.csv — rank per sample, mean WT/G3 rank, rank_delta.
Task 3: CFS overlap for every gene in the top-30 union.

Run from the repo root:
    python3 scripts/u2os_annotate_rankshift_cfs.py
"""

import glob
import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR   = Path("data/u2os")
CFS_DIR    = Path("data/CFS_fixed")
RESULT_DIR = Path("results/u2os")

SAMPLES    = ["WT-U", "WT-A", "G3-U", "G3-A"]

# GRCh38 PAR coordinates (both X and Y copies give the same gene,
# so we flag by gene name derived from coordinate check below)
PAR1_X = (10_001, 2_781_479)
PAR2_X = (155_701_383, 156_030_895)
PAR1_Y = (10_001, 2_781_479)
PAR2_Y = (56_887_903, 57_217_415)


# ── Step 0: extract gene coordinates from the TSV files ───────────────────────
# We pull CHROM + POS for every PICK1 protein-coding variant and take
# min/max POS per gene as approximate gene span. We only need the union
# of top-30 genes, so we filter early.

print("=== Step 0: Extracting gene coordinates from TSVs ===\n")

# Build target gene set first (cheap — just reads existing CSVs)
target_genes: set[str] = set()
top30: dict[str, pd.DataFrame] = {}
for s in SAMPLES:
    df = pd.read_csv(RESULT_DIR / f"top30_{s}.csv")
    top30[s] = df
    target_genes |= set(df["Gene"])

print(f"Union of top-30 genes: {len(target_genes)}")

# One-pass scan per TSV — keep only target genes
gene_coords: dict[str, dict] = {}   # gene -> {chrom, start, end}

for s in SAMPLES:
    path = DATA_DIR / f"{s}_vep_anno_non_common.tsv.gz"
    print(f"Scanning {s} for coordinates ...", flush=True)
    reader = pd.read_csv(
        path, sep="\t",
        usecols=["CHROM", "POS", "PICK", "BIOTYPE", "SYMBOL"],
        low_memory=False,
        chunksize=200_000,
    )
    for chunk in reader:
        sub = chunk[
            (chunk["PICK"] == 1) &
            (chunk["BIOTYPE"] == "protein_coding") &
            (chunk["SYMBOL"].isin(target_genes))
        ]
        for gene, grp in sub.groupby("SYMBOL"):
            chrom = grp["CHROM"].iloc[0]
            mn    = int(grp["POS"].min())
            mx    = int(grp["POS"].max())
            if gene not in gene_coords:
                gene_coords[gene] = {"chrom": chrom, "start": mn, "end": mx}
            else:
                gene_coords[gene]["start"] = min(gene_coords[gene]["start"], mn)
                gene_coords[gene]["end"]   = max(gene_coords[gene]["end"],   mx)

print(f"\nCoordinates resolved for {len(gene_coords)} / {len(target_genes)} genes\n")

# Genes with no PICK1+protein_coding hits (should not happen, but guard)
missing = target_genes - set(gene_coords)
if missing:
    print(f"WARNING — no coordinates found for: {missing}\n")


# ── Step 1: flag derivation ────────────────────────────────────────────────────
print("=== Step 1: Flagging PAR artefacts and HLA locus genes ===\n")

def classify_flag(gene: str) -> str:
    if gene.startswith("HLA-"):
        return "HLA_locus"
    if gene not in gene_coords:
        return ""
    chrom = gene_coords[gene]["chrom"]
    pos   = (gene_coords[gene]["start"] + gene_coords[gene]["end"]) // 2
    if chrom == "chrX":
        if PAR1_X[0] <= pos <= PAR1_X[1]:
            return "PAR_artefact"
        if PAR2_X[0] <= pos <= PAR2_X[1]:
            return "PAR_artefact"
    if chrom == "chrY":
        if PAR1_Y[0] <= pos <= PAR1_Y[1]:
            return "PAR_artefact"
        if PAR2_Y[0] <= pos <= PAR2_Y[1]:
            return "PAR_artefact"
    return ""

flag_map = {gene: classify_flag(gene) for gene in target_genes}

flagged = [(g, f) for g, f in flag_map.items() if f]
print("Flagged genes:")
for g, f in sorted(flagged):
    c = gene_coords.get(g, {})
    print(f"  {g:15s}  {f:15s}  {c.get('chrom','?')}:{c.get('start','?')}-{c.get('end','?')}")
print()

# Write updated CSVs
for s in SAMPLES:
    df = top30[s].copy()
    df["flag"] = df["Gene"].map(flag_map).fillna("")
    out = RESULT_DIR / f"top30_{s}.csv"
    df.to_csv(out, index=False)
    print(f"  Updated → {out}")

print()


# ── Step 2: rank-shift table ───────────────────────────────────────────────────
print("=== Step 2: Building rank_shift_u2os.csv ===\n")

# rank lookup: gene -> sample -> rank (1-based); absent = NaN
rank_lookup: dict[str, dict[str, float]] = {g: {} for g in target_genes}
for s in SAMPLES:
    for _, row in top30[s].iterrows():
        rank_lookup[row["Gene"]][s] = float(row["Rank"])

rows = []
for gene in sorted(target_genes):
    r = {s: rank_lookup[gene].get(s, np.nan) for s in SAMPLES}

    # For averaging, treat absent (>30) as 31
    wt_vals = [r["WT-U"] if not np.isnan(r["WT-U"]) else 31.0,
               r["WT-A"] if not np.isnan(r["WT-A"]) else 31.0]
    g3_vals = [r["G3-U"] if not np.isnan(r["G3-U"]) else 31.0,
               r["G3-A"] if not np.isnan(r["G3-A"]) else 31.0]

    mean_wt = np.mean(wt_vals)
    mean_g3 = np.mean(g3_vals)
    delta   = mean_wt - mean_g3   # positive = rises in G3 (enriched in KO)

    rows.append({
        "Gene":          gene,
        "flag":          flag_map[gene],
        "WT-U_rank":     int(r["WT-U"]) if not np.isnan(r["WT-U"]) else ">30",
        "WT-A_rank":     int(r["WT-A"]) if not np.isnan(r["WT-A"]) else ">30",
        "G3-U_rank":     int(r["G3-U"]) if not np.isnan(r["G3-U"]) else ">30",
        "G3-A_rank":     int(r["G3-A"]) if not np.isnan(r["G3-A"]) else ">30",
        "mean_WT_rank":  round(mean_wt, 1),
        "mean_G3_rank":  round(mean_g3, 1),
        "rank_delta":    round(delta,   1),
    })

rank_df = pd.DataFrame(rows).sort_values("rank_delta", ascending=False)
out_rank = RESULT_DIR / "rank_shift_u2os.csv"
rank_df.to_csv(out_rank, index=False)
print(f"Saved → {out_rank}\n")

print("Top 15 by rank_delta (most enriched in G3 KO):")
print(rank_df.head(15).to_string(index=False))
print()
print("Bottom 15 by rank_delta (most enriched in WT):")
print(rank_df.tail(15).to_string(index=False))
print()


# ── Step 3: CFS overlap ────────────────────────────────────────────────────────
print("=== Step 3: CFS overlap for top-30 union ===\n")

# Load all CFS BED files
bed_files = list(CFS_DIR.glob("*.bed"))
if not bed_files:
    raise FileNotFoundError(f"No BED files found in {CFS_DIR}")

cfs_df = pd.concat(
    [pd.read_csv(b, sep="\t", header=None,
                 names=["chr", "start", "end", "id", "score", "strand"])
     for b in bed_files],
    ignore_index=True,
)
print(f"Loaded {len(cfs_df)} CFS regions from {len(bed_files)} BED files")

def cfs_hits(gene: str) -> list[str]:
    if gene not in gene_coords:
        return []
    c   = gene_coords[gene]
    sub = cfs_df[cfs_df["chr"] == c["chrom"]]
    ov  = sub[(c["start"] <= sub["end"]) & (c["end"] >= sub["start"])]
    return ov["id"].tolist()

cfs_rows = []
for s in SAMPLES:
    for _, row in top30[s].iterrows():
        gene  = row["Gene"]
        hits  = cfs_hits(gene)
        cfs_rows.append({
            "Sample":      s,
            "Rank":        int(row["Rank"]),
            "Gene":        gene,
            "flag":        flag_map[gene],
            "n_cfs":       len(hits),
            "CFS_regions": ",".join(hits) if hits else "",
        })

cfs_out = pd.DataFrame(cfs_rows)
out_cfs = RESULT_DIR / "cfs_overlap_top30_u2os.csv"
cfs_out.to_csv(out_cfs, index=False)
print(f"Saved → {out_cfs}\n")

# Summary per sample
print("── CFS overlap summary per sample ──")
for s in SAMPLES:
    sub   = cfs_out[cfs_out["Sample"] == s]
    hits  = sub[sub["n_cfs"] > 0]
    print(f"\n{s}: {len(hits)}/30 top genes overlap a CFS")
    for _, r in hits.iterrows():
        flag_str = f"  [{r['flag']}]" if r["flag"] else ""
        print(f"  rank {r['Rank']:2d}  {r['Gene']:15s}  {r['CFS_regions']}{flag_str}")

# Spotlight genes
print("\n── Spotlight genes ──")
for spotlight in ["FHIT", "ROBO2"]:
    rows_sp = cfs_out[cfs_out["Gene"] == spotlight]
    if rows_sp.empty:
        print(f"  {spotlight}: not in any top-30")
    else:
        for _, r in rows_sp.iterrows():
            status = r["CFS_regions"] if r["CFS_regions"] else "no CFS overlap"
            print(f"  {spotlight}  {r['Sample']}  rank {r['Rank']:2d}  {status}")

print("\nDone.")
