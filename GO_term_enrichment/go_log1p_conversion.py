#!/usr/bin/env python3
"""
Convert a GO-term x genome count table to log-transformed values for heatmap visualization.

INPUT format (like go_counts_Pilosi.csv):
- Column 1: GO term / feature name
- Columns 2..N: genome names
- Values: raw gene counts (integers)

Transformation:
- Natural log1p: ln(1 + x)  (works with zeros)

USAGE:
  python go_log1p_table.py input.csv
  python go_log1p_table.py input.xlsx

OUTPUT:
  log_out/
    <input_basename>_log1p.csv
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    p = path.lower()
    if p.endswith(".csv"):
        return pd.read_csv(path)
    if p.endswith(".xlsx") or p.endswith(".xls"):
        return pd.read_excel(path)
    raise ValueError("Input must be .csv, .xlsx, or .xls")


def main():
    if len(sys.argv) != 2:
        print("USAGE: python go_log1p_table.py input.(csv|xlsx|xls)")
        sys.exit(1)

    in_path = sys.argv[1]
    df = read_table(in_path)

    if df.shape[1] < 3:
        raise ValueError("Table must have >= 1 GO-term column + >= 2 genome columns.")

    term_col = df.columns[0]
    genome_cols = [c for c in df.columns if c != term_col]

    # Clean: treat empty strings as missing, fill missing with 0
    df[genome_cols] = df[genome_cols].replace("", np.nan).fillna(0)

    # Ensure numeric and non-negative
    vals = df[genome_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    if (vals < 0).any().any():
        raise ValueError("Found negative values. Raw counts should be >= 0.")

    # log1p transform
    df_log = df.copy()
    df_log[genome_cols] = np.log1p(vals)

    outdir = Path("log_out")
    outdir.mkdir(exist_ok=True)

    base = Path(in_path).stem
    out_file = outdir / f"{base}_log1p.csv"
    df_log.to_csv(out_file, index=False)

    print(f"[OK] log1p table written to: {out_file.resolve()}")
    print(f"Feature column: {term_col}")
    print(f"Genomes (n={len(genome_cols)}): {', '.join(map(str, genome_cols[:8]))}"
          + (" ..." if len(genome_cols) > 8 else ""))


if __name__ == "__main__":
    main()
