#!/usr/bin/env python3
import os
import re
import gc
import glob
import argparse
import pandas as pd

COLS = [

    "chrom", "start", "end", "modification", "score", "strand",
    "start_c", "end_c", "color", "Nvalid_cov", "fraction_modified",
    "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"
]

def main():
    parser = argparse.ArgumentParser(description="Preprocess modkit bedMethyl files")
    parser.add_argument("--input-dir", required=True, help="Directory containing .bed files")
    parser.add_argument("--output-dir", required=True, help="Directory for outputs")
    parser.add_argument("--min-coverage", type=int, default=5, help="Minimum Nvalid_cov threshold")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    min_coverage = args.min_coverage
    per_sample_dir = os.path.join(output_dir, "per_sample_filtered")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(per_sample_dir, exist_ok=True)

    bed_files = sorted(glob.glob(os.path.join(input_dir, "*.bed")))
    if not bed_files:
        raise FileNotFoundError(f"No .bed files found in: {input_dir}")

    print(f"Found {len(bed_files)} bed files")

    for f in bed_files:
        print(" -", os.path.basename(f))

    summary_rows = []
    filtered_paths = []

    # process each sample separately first

    for bed_path in bed_files:
        sample = os.path.splitext(os.path.basename(bed_path))[0]
        age_match = re.search(r"(\d+)", sample)
        age = int(age_match.group(1)) if age_match else pd.NA

        print(f"\nProcessing {sample} ...")
        df = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            names=COLS,
            low_memory=False
        )

        raw_rows = len(df)

        df["sample"] = sample
        df["age"] = age
        df["ont_cpg"] = df["chrom"].astype(str) + "_" + df["start"].astype(str)

        df = df[df["Nvalid_cov"] >= min_coverage].copy()
        filtered_rows = len(df)

        out_file = os.path.join(per_sample_dir, f"{sample}_filtered.tsv")
        df.to_csv(out_file, sep="\t", index=False)
        filtered_paths.append(out_file)

        summary_rows.append({
            "sample": sample,
            "age": age,
            "raw_rows": raw_rows,
            "filtered_rows": filtered_rows
        })

        print(f"{sample}: raw={raw_rows:,}, filtered={filtered_rows:,}")
        del df
        gc.collect()

    # save processing summary
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(output_dir, "filtering_summary.tsv"), sep="\t", index=False)

    print("\nConcatenating filtered files ...")
    all_data = []
    for path in filtered_paths:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        all_data.append(df)

    pbr = pd.concat(all_data, ignore_index=True)
    print(f"Combined pbr shape: {pbr.shape}")

    pbr_out = os.path.join(output_dir, "pbr_filtered.tsv")
    pbr.to_csv(pbr_out, sep="\t", index=False)

    # split by modification type
    ont_h = pbr[pbr["modification"] == "h"].copy()
    ont_m = pbr[pbr["modification"] == "m"].copy()
    hm_keys = (
        pbr.groupby("ont_cpg")["modification"]
        .apply(lambda x: "".join(sorted(set(x))))
        .reset_index(name="types")
    )

    hm_keys = hm_keys[hm_keys["types"] == "hm"]["ont_cpg"]
    ont_hm = pbr[pbr["ont_cpg"].isin(hm_keys)].copy()

    print(f"ont_h shape:  {ont_h.shape}")
    print(f"ont_m shape:  {ont_m.shape}")
    print(f"ont_hm shape: {ont_hm.shape}")

    ont_h.to_csv(os.path.join(output_dir, "ont_h.tsv"), sep="\t", index=False)
    ont_m.to_csv(os.path.join(output_dir, "ont_m.tsv"), sep="\t", index=False)
    ont_hm.to_csv(os.path.join(output_dir, "ont_hm.tsv"), sep="\t", index=False)

    print("\nDone.")
    print(f"Outputs written to: {output_dir}")

if __name__ == "__main__":
    main()
