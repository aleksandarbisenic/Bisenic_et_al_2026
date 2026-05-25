#!/usr/bin/env python3
# differential_go_terms.py
#
# PURPOSE
#   Count (per genome) how many UNIQUE genes are annotated with each GO description
#   listed in go_descriptions.txt, and write go_counts_output.csv.
#
# INPUT MODES
#   1) Excel/CSV wide paired columns (recommended; same as Fisher_GO_input.xlsx):
#        <Genome1 gene_id> | GO term | <Genome2 gene_id> | GO term.1 | ...
#      Example:
#        python differential_go_terms.py go_descriptions.txt Fisher_GO_input.xlsx
#
#   2) Legacy per-genome TXT/TSV files (backwards compatible with the old script):
#      Assumes the 4th column (index 3) contains GO descriptions.
#        python differential_go_terms.py go_descriptions.txt genome1.txt genome2.txt ...
#
# OUTPUT
#   go_counts_output.csv with columns:
#     GO Description, <genome1>, <genome2>, ...
#
# NOTES
#   - Matching uses substring containment (same spirit as the original script),
#     but with regex disabled to avoid surprises from special characters.
#   - Counting is by UNIQUE gene_id per genome for each GO description. This matches the
#     logic used in go_enrichment_fisher.py when building counts.

import argparse
import os
import sys
import re
from typing import Dict, List, Tuple

import pandas as pd


def _normalize_go_desc(s: str) -> str:
    """Normalize GO description strings for exact matching (trim + collapse whitespace)."""
    if s is None:
        return ""
    s = str(s).strip()
    s = re.sub(r"\s+", " ", s)
    return s



def read_go_descriptions(file_path: str) -> List[str]:
    """Read GO descriptions (one per line)."""
    with open(file_path, "r", encoding="utf-8") as f:
        return [_normalize_go_desc(line) for line in f if line.strip()]


def _load_table(path: str, sheet: str | None = None) -> pd.DataFrame:
    pl = path.lower()
    if pl.endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(path, sheet_name=sheet or 0)
    if pl.endswith(".csv"):
        return pd.read_csv(path)
    if pl.endswith(".tsv"):
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported input type for {path}. Use .xlsx/.xls/.csv/.tsv or legacy per-genome .txt/.tsv.")


def to_long_from_wide(df: pd.DataFrame) -> pd.DataFrame:
    """Convert wide paired columns into long format: genome, gene_id, go_term."""
    cols = list(df.columns)
    rows = []
    i = 0

    # Pair strategy: <gene col>, <GO col> where GO col name starts with 'go term' (case-insensitive)
    while i < len(cols) - 1:
        gene_col = cols[i]
        go_col = cols[i + 1]
        if str(go_col).strip().lower().startswith("go term"):
            genome = str(gene_col).strip()
            sub = df[[gene_col, go_col]].dropna(subset=[gene_col, go_col]).copy()
            sub.columns = ["gene_id", "go_term"]
            sub["genome"] = genome
            rows.append(sub[["genome", "gene_id", "go_term"]])
            i += 2
        else:
            i += 1

    if not rows:
        raise ValueError(
            "Could not parse paired columns. Expected alternating '<genome gene_id column>', 'GO term' columns."
        )

    long_df = pd.concat(rows, ignore_index=True)

    # Standardize types and drop exact duplicates
    long_df["genome"] = long_df["genome"].astype(str)
    long_df["gene_id"] = long_df["gene_id"].astype(str)
    long_df["go_term"] = long_df["go_term"].astype(str)
    long_df = long_df.dropna(subset=["genome", "gene_id", "go_term"]).drop_duplicates()

    return long_df


def count_from_long(go_descriptions: List[str], long_df: pd.DataFrame) -> pd.DataFrame:
    """Return a dataframe: GO Description x genomes with unique-gene counts."""
    # Ensure required columns exist
    for c in ("genome", "gene_id", "go_term"):
        if c not in long_df.columns:
            raise ValueError(f"Long dataframe missing required column: {c}")

    # Deduplicate gene-go assignments
    long_df = long_df.drop_duplicates(subset=["genome", "gene_id", "go_term"]).copy()

    # Prepare output
    genomes = sorted(long_df["genome"].unique().tolist())
    out = pd.DataFrame({"GO Description": go_descriptions})

    # For speed: pre-split by genome
    by_genome = {g: long_df.loc[long_df["genome"] == g, ["gene_id", "go_term"]].copy() for g in genomes}

    for g in genomes:
        sub = by_genome[g]
        # For each description, count unique genes with an EXACT GO-term description match
        counts = []
        go_series = sub["go_term"].astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
        gene_series = sub["gene_id"].astype(str)

        for desc in go_descriptions:
            # Exact match after normalization (trim + collapse whitespace)
            mask = go_series.eq(desc)
            counts.append(int(gene_series[mask].nunique()))
        out[g] = counts

    return out


def count_legacy_file(go_descriptions: List[str], bacterial_file: str) -> Dict[str, int]:
    """Legacy mode: bacterial_file is a TSV-like file; 4th column (index 3) contains GO descriptions."""
    df = pd.read_csv(bacterial_file, sep="\t", header=None, dtype=str)
    if df.shape[1] <= 3:
        raise ValueError(f"Legacy file {bacterial_file} has <4 columns; expected GO descriptions in column 4 (index 3).")

    go_terms = df.iloc[:, 3].astype(str)
    # Old script counted rows; keep that behavior for legacy files
    go_terms_norm = go_terms.astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
    return {desc: int((go_terms_norm == desc).sum()) for desc in go_descriptions}


def main():
    ap = argparse.ArgumentParser(
        description="Count per-genome GO-description gene counts (supports Fisher_GO_input.xlsx wide format)."
    )
    ap.add_argument("go_descriptions", help="Path to go_descriptions.txt (one GO description per line)")
    ap.add_argument(
        "inputs",
        nargs="+",
        help="Either a single wide Excel/CSV (e.g., Fisher_GO_input.xlsx) OR legacy per-genome TXT/TSV files",
    )
    ap.add_argument("--sheet", default=None, help="Excel sheet name (if Excel input). Default: first sheet.")
    ap.add_argument("--out", default="go_counts_output.csv", help="Output CSV name (default: go_counts_output.csv)")
    args = ap.parse_args()

    go_descriptions = read_go_descriptions(args.go_descriptions)
    if not go_descriptions:
        sys.exit("No GO descriptions found (file was empty after stripping blank lines).")

    # If they provided exactly one spreadsheet-like file => wide mode
    if len(args.inputs) == 1 and args.inputs[0].lower().endswith((".xlsx", ".xlsm", ".xls", ".csv", ".tsv")):
        df = _load_table(args.inputs[0], sheet=args.sheet)
        long_df = to_long_from_wide(df)
        out_df = count_from_long(go_descriptions, long_df)
        out_df.to_csv(args.out, index=False)
        print(f"[OK] Wrote {args.out} (input: {args.inputs[0]}; genomes: {out_df.shape[1]-1}).")
        return

    # Otherwise: legacy per-genome files
    output_df = pd.DataFrame({"GO Description": go_descriptions})
    for bacterial_file in args.inputs:
        col_name = os.path.basename(bacterial_file)
        counts = count_legacy_file(go_descriptions, bacterial_file)
        output_df[col_name] = output_df["GO Description"].map(counts).astype(int)

    output_df.to_csv(args.out, index=False)
    print(f"[OK] Wrote {args.out} (legacy mode; files: {len(args.inputs)}).")


if __name__ == "__main__":
    main()
