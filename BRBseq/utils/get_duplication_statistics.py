import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


def load_matrix(raw_dir, filename):
    return mmread(Path(raw_dir) / filename).tocsc()  # Convert to CSC format for efficient column operations


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate UMI deduplication statistics.")
    parser.add_argument("--raw_dir", type=str, required=True, help="STAR Directory containing raw matrices and barcodes.")
    parser.add_argument(
        "--dedup_method",
        type=str,
        required=True,
        help="Deduplication method used (e.g., '1MM_Directional').",
    )
    parser.add_argument("--output_csv", type=str, required=True, help="Path to output CSV file for duplication statistics.")
    parser.add_argument(
        "--manifest",
        type=str,
        default=None,
        help="Optional two-column TSV (sample_id, barcode) used to annotate barcodes with sample IDs.",
    )
    args = parser.parse_args()

    # Load barcode and feature names
    barcodes = pd.read_csv(Path(args.raw_dir) / "barcodes.tsv", header=None)[0].tolist()

    # Load matrices
    mat_umi = load_matrix(args.raw_dir, f"umiDedup-{args.dedup_method}.mtx")  # deduplicated
    mat_reads = load_matrix(args.raw_dir, "umiDedup-NoDedup.mtx")  # raw reads

    if mat_umi.shape != mat_reads.shape:
        sys.exit("ERROR: UMI and NoDedup matrices have different shapes.")

    # colSums — equivalent to R colSums(), summing over all genes per barcode
    umi_counts = np.asarray(mat_umi.sum(axis=0)).flatten()  # shape: (n_barcodes,)
    read_counts = np.asarray(mat_reads.sum(axis=0)).flatten()

    with np.errstate(divide="ignore", invalid="ignore"):
        umi_pct = np.where(read_counts > 0, np.round(umi_counts / read_counts * 100, 2), np.nan)

    result = pd.DataFrame({
        "barcode": barcodes,
        "umi_count": umi_counts.astype(int),
        "read_count": read_counts.astype(int),
        "umi_pct": umi_pct,
    })

    if args.manifest is not None:
        wl = pd.read_csv(args.manifest, sep="\t", header=0, names=["sample_id", "barcode"])
        result = result.merge(wl[["barcode", "sample_id"]], on="barcode", how="left")
        result.insert(0, "sample_id", result.pop("sample_id"))
        result.sort_values(by="sample_id", inplace=True)

    summary = pd.DataFrame({
        "total_barcodes": [len(barcodes)],
        "total_umi": [umi_counts.sum()],
        "total_reads": [read_counts.sum()],
        "overall_umi_pct": [
            np.round(umi_counts.sum() / read_counts.sum() * 100, 2) if read_counts.sum() > 0 else np.nan
        ],
    })

    result.to_csv(Path(args.output_csv), index=False)
    summary.to_csv(Path(args.output_csv).with_suffix(".summary.csv"), index=False)