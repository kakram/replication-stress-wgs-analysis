import pandas as pd
import zipfile
import glob
import os

# --------------------------------------------------
# 1. Define directories
# --------------------------------------------------

BASE_DIR = os.path.dirname(os.path.abspath(__file__))  # /scripts
PROJECT_ROOT = os.path.dirname(BASE_DIR)               # directory above /scripts

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
CFS_EXTRACT_DIR = os.path.join(RESULTS_DIR, "cfs_extracted")

# Make sure output directories exist
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(CFS_EXTRACT_DIR, exist_ok=True)

print("Data directory:", DATA_DIR)
print("Results directory:", RESULTS_DIR)


# --------------------------------------------------
# 2. Load mutated genes file
# --------------------------------------------------

genes_file = os.path.join(DATA_DIR, "MCF7MutatedGenesLocations.xlsx")
genes_raw = pd.read_excel(genes_file, header=1)

# Remove placeholder column if present
if "Unnamed: 0" in genes_raw.columns:
    genes_raw = genes_raw.drop(columns=["Unnamed: 0"])

genes_raw["Gene_Start"] = pd.to_numeric(genes_raw["Gene_Start"])
genes_raw["Gene_End"] = pd.to_numeric(genes_raw["Gene_End"])
genes_raw["CHROMOSOME"] = genes_raw["CHROMOSOME"].astype(str)

print("Loaded genes:", genes_raw.shape)


# --------------------------------------------------
# 3. Extract CFS BED files from ZIP
# --------------------------------------------------

cfs_zip_file = os.path.join(DATA_DIR, "cfsfixed.zip")

with zipfile.ZipFile(cfs_zip_file, "r") as z:
    z.extractall(CFS_EXTRACT_DIR)

print("Extracted BED files to:", CFS_EXTRACT_DIR)


# --------------------------------------------------
# 4. Load all extracted BED files
# --------------------------------------------------

bed_files = glob.glob(os.path.join(CFS_EXTRACT_DIR, "*.bed"))

cfs_list = []
for bf in bed_files:
    df = pd.read_csv(
        bf, sep="\t", header=None,
        names=["chr", "start", "end", "id", "score", "strand"]
    )
    cfs_list.append(df)

cfs = pd.concat(cfs_list, ignore_index=True)

print("Loaded CFS entries:", cfs.shape)


# --------------------------------------------------
# 5. CFS overlap computation
# --------------------------------------------------

def find_cfs_overlap(row):
    """Return comma-separated fragile site IDs overlapped by this gene."""
    chr_name = row["CHROMOSOME"]
    sub = cfs[cfs["chr"] == chr_name]

    overlaps = sub[
        (row["Gene_Start"] <= sub["end"]) &
        (row["Gene_End"] >= sub["start"])
    ]

    if overlaps.empty:
        return ""

    return ",".join(overlaps["id"].tolist())


genes_raw["CFS_overlap"] = genes_raw.apply(find_cfs_overlap, axis=1)

print("Overlap computation complete.")


# --------------------------------------------------
# 6. Save output file
# --------------------------------------------------

output_file = os.path.join(RESULTS_DIR, "MutatedGenes_CFS_annotation_fresh.xlsx")
genes_raw.to_excel(output_file, index=False)

print("Saved output to:", output_file)


