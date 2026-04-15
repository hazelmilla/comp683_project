import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

# ============================================================
# 1. LOAD DATA
# ============================================================

print("Loading data...")

df = pd.read_csv(
    "/work/users/d/p/dpguilba/meth_clock/unfiltered/ont_m.tsv",
    sep="\t",
    dtype={
        "chrom": "category",
        "sample": "category",
        "ont_cpg": "category"
    }
)

print("Shape:", df.shape)

# ============================================================
# 2. BASIC FEATURES
# ============================================================

df["coverage"] = df["Nmod"] + df["Ncanonical"]
df = df[df["coverage"] > 0]

df["beta"] = df["Nmod"] / df["coverage"]

# ============================================================
# 3. GLOBAL CpG SET (CRITICAL)
# ============================================================

print("Building global CpG set...")

global_cpgs = df[["chrom", "start", "ont_cpg"]].drop_duplicates()
global_cpgs = global_cpgs.sort_values(["chrom", "start"]).reset_index(drop=True)

print("Total unique CpGs:", len(global_cpgs))

# ============================================================
# 4. TRUE IMPUTATION FUNCTION
# ============================================================

def process_sample_full(sample_df, global_cpgs, frac=0.01):

    # merge onto full CpG set
    sample_df = global_cpgs.merge(
        sample_df[["chrom", "start", "ont_cpg", "beta"]],
        on=["chrom", "start", "ont_cpg"],
        how="left"
    )

    out = []

    for chrom, chrom_df in sample_df.groupby("chrom"):

        chrom_df = chrom_df.sort_values("start")

        observed = chrom_df.dropna(subset=["beta"])

        # skip chromosomes with too little data
        if len(observed) < 50:
            chrom_df["beta_loess"] = np.nan
            chrom_df["beta_imputed"] = chrom_df["beta"]
            out.append(chrom_df)
            continue

        # fit LOESS on observed CpGs only
        smoothed = lowess(
            observed["beta"],
            observed["start"],
            frac=frac,
            it=0,
            return_sorted=True
        )

        x_smooth = smoothed[:, 0]
        y_smooth = smoothed[:, 1]

        # interpolate to ALL CpGs (including missing)
        beta_loess_full = np.interp(
            chrom_df["start"],
            x_smooth,
            y_smooth
        )

        chrom_df = chrom_df.copy()
        chrom_df["beta_loess"] = beta_loess_full

        # impute ONLY missing sites
        chrom_df["beta_imputed"] = chrom_df["beta"]

        missing = chrom_df["beta"].isna()
        chrom_df.loc[missing, "beta_imputed"] = chrom_df.loc[missing, "beta_loess"]

        out.append(chrom_df)

    return pd.concat(out, ignore_index=True)

# ============================================================
# 5. PROCESS EACH SAMPLE (memory safe)
# ============================================================

print("Imputing each sample...")

wide_data = []   # will store series for each sample

samples = df["sample"].unique()

for sample in samples:
    print("Processing:", sample)

    sdf = df[df["sample"] == sample]

    processed = process_sample_full(sdf, global_cpgs, frac=0.01)

    # keep only needed columns for wide matrix
    series = processed.set_index("ont_cpg")["beta_imputed"]
    series.name = sample

    wide_data.append(series)

# ============================================================
# 6. BUILD WIDE MATRIX
# ============================================================

print("Building wide matrix...")

wide_df = pd.concat(wide_data, axis=1)

print("Wide shape:", wide_df.shape)

# ============================================================
# 7. SAVE OUTPUT
# ============================================================

# Save efficient version
wide_df.to_parquet("ont_m_imputed_wide.parquet")

# Optional CSV (can be very large!)
# wide_df.to_csv("ont_m_imputed_wide.csv")

print("Done.")