import pandas as pd
import zipfile
import glob
import os

# --------------------------------------------------
# Directories
# --------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))     # /scripts
ROOT_DIR = os.path.dirname(BASE_DIR)                      # project root

DATA_DIR = os.path.join(ROOT_DIR, "data")
RESULTS_DIR = os.path.join(ROOT_DIR, "results")
CFS_DIR = os.path.join(RESULTS_DIR, "cfs_extracted")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(CFS_DIR, exist_ok=True)

# --------------------------------------------------
# Load mutated genes
# --------------------------------------------------
genes_path = os.path.join(DATA_DIR, "MCF7MutatedGenesLocations.xlsx")
genes = pd.read_excel(genes_path, header=1)

# Remove placeholder column if present
if "Unnamed: 0" in genes.columns:
    genes = genes.drop(columns=["Unnamed: 0"])

genes["Gene_Start"] = pd.to_numeric(genes["Gene_Start"])
genes["Gene_End"] = pd.to_numeric(genes["Gene_End"])
genes["CHROMOSOME"] = genes["CHROMOSOME"].astype(str)

# --------------------------------------------------
# Load CFS regions (BED files inside ZIP)
# --------------------------------------------------
cfs_zip = os.path.join(DATA_DIR, "cfsfixed.zip")
with zipfile.ZipFile(cfs_zip, "r") as z:
    z.extractall(CFS_DIR)

cfs_list = []
for bed in glob.glob(os.path.join(CFS_DIR, "*.bed")):
    df = pd.read_csv(
        bed,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "id", "score", "strand"]
    )
    cfs_list.append(df)

cfs = pd.concat(cfs_list, ignore_index=True)

# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def overlaps_cfs(row):
    """True if gene intersects any fragile site."""
    sub = cfs[cfs["chr"] == row["CHROMOSOME"]]
    ov = sub[(row["Gene_Start"] <= sub["end"]) &
             (row["Gene_End"] >= sub["start"])]
    return not ov.empty

def nearest_cfs_distance(row):
    """Return distance and CFS ID for nearest fragile site. 
       Returns (None, None) if chromosome has no CFS."""
    sub = cfs[cfs["chr"] == row["CHROMOSOME"]]

    if sub.empty:
        # No fragile sites on this chromosome
        return (None, None)

    distances = []

    for _, cs in sub.iterrows():
        if row["Gene_End"] < cs["start"]:
            d = cs["start"] - row["Gene_End"]
        elif row["Gene_Start"] > cs["end"]:
            d = row["Gene_Start"] - cs["end"]
        else:
            d = 0
        distances.append((d, cs["id"]))

    # Now guaranteed not empty
    return min(distances, key=lambda x: x[0])

# --------------------------------------------------
# Compute distances for non-overlapping genes
# --------------------------------------------------
genes["overlaps_CFS"] = genes.apply(overlaps_cfs, axis=1)

non_overlapping = genes[~genes["overlaps_CFS"]].copy()

non_overlapping[["distance_from_CFS", "nearest_CFS"]] = (
    non_overlapping.apply(nearest_cfs_distance, axis=1, result_type="expand")
)

# --------------------------------------------------
# Extract genes within 1 Mb (CFS-adjacent)
# --------------------------------------------------
NEAR_THRESHOLD = 1_000_000  # 1 Mb

near_genes = non_overlapping[
    non_overlapping["distance_from_CFS"] < NEAR_THRESHOLD
].copy()

# --------------------------------------------------
# Save output table
# --------------------------------------------------
output_path = os.path.join(RESULTS_DIR, "Genes_near_CFS_within_1Mb.xlsx")
near_genes.to_excel(output_path, index=False)

print("Saved:", output_path)
print("Number of genes within 1 Mb of a CFS:", near_genes.shape[0])

