# Script to build feature matrix for kNN/Leiden graphing using drift-associated CpGs and imputed methylation data
#!/usr/bin/env python3
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd

BASE = Path("/work/users/k/a/katlyssa/Sample BED files")
PARQUET = BASE / "full_impute.parquet"
DRIFT = BASE / "driftCpG_robust_coef.csv"
OUT = BASE / "drift_feature_matrix.tsv"
WHITE_P = 0.05
TOP_N = 2000   # adjust: 500 / 1000 / 2000

print("Loading drift data...")

drift = pd.read_csv(DRIFT)

# Filter drift CpGs
drift = drift[drift["whiteP"] < WHITE_P].copy()

# Fix CpG formatting (remove .0 issue)
drift["ont_cpg"] = drift["ont_cpg"].astype(str)

# Rank by drift magnitude
drift["abs_drift"] = drift["drift_std"].abs()

drift = drift.sort_values("abs_drift", ascending=False)

if TOP_N:
    drift = drift.head(TOP_N)

drift_set = set(drift["ont_cpg"])

print(f"Drift CpGs retained: {len(drift_set)}")

print("Reading parquet in chunks...")

pf = pq.ParquetFile(PARQUET)

dfs = []

for i in range(pf.num_row_groups):

    print(f"Reading row group {i+1}/{pf.num_row_groups}")

    table = pf.read_row_group(i)
    df = table.to_pandas()

    # Fix CpG naming
    df["ont_cpg"] = df["ont_cpg"].astype(str).str.replace(".0", "", regex=False)

    # Filter only drift CpGs
    df = df[df["ont_cpg"].isin(drift_set)]
  
    if len(df) > 0:
        dfs.append(df)

print("Concatenating filtered data...")

df_all = pd.concat(dfs, ignore_index=True)

print("Final shape:", df_all.shape)

# Select methylation columns
meth_cols = [c for c in df_all.columns if c.startswith("meth_basic_")]

print("Methylation columns:")
print(meth_cols)

# Build matrix
mat = df_all[["ont_cpg"] + meth_cols].copy()

# Remove duplicates (important)
mat = mat.drop_duplicates(subset="ont_cpg")

mat = mat.set_index("ont_cpg")

print("Matrix shape (CpGs x samples):", mat.shape)

# Save
mat.to_csv(OUT, sep="\t")



print(f"Saved matrix to: {OUT}")
