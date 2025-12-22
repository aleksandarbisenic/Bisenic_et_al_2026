#!/usr/bin/env python3
# differential_go_terms_exact.py
# Usage:
#   python differential_go_terms_exact.py go_descriptions.txt file1.tsv file2.tsv ...
#
# Common options:
#   --sep "\t"               # input sep (default: tab)
#   --header -1              # no header (default). Use 0 if your file has a header row
#   --col 3                  # zero-based index of GO column (default: 3 => 4th column)
#   --split-regex "\s*[;|,]\s*"   # split multiple terms per cell; compare exact tokens
#   --ignore-case            # make matching case-insensitive (default: case-sensitive)
#   --collapse-spaces        # collapse internal whitespace in terms (default: on)
#   --out go_counts_output.csv

import argparse
import os
import sys
import re
import pandas as pd

def read_go_descriptions(path, ignore_case=False, collapse_spaces=True):
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f.readlines()]
    # normalize
    if collapse_spaces:
        lines = [re.sub(r"\s+", " ", x) for x in lines]
    if ignore_case:
        lines = [x.lower() for x in lines]
    # keep order, but also have a set for quick lookup
    return lines

def load_go_column(bacterial_file, sep, header, col_idx):
    df = pd.read_csv(bacterial_file, sep=sep, header=None if header == -1 else header, dtype=str)
    try:
        series = df.iloc[:, col_idx].astype(str)
    except Exception:
        raise IndexError(f"Column index {col_idx} out of range for file: {bacterial_file}")
    return series

def normalize_series(ser: pd.Series, ignore_case=False, collapse_spaces=True):
    ser = ser.fillna("").astype(str).str.strip()
    if collapse_spaces:
        ser = ser.str.replace(r"\s+", " ", regex=True)
    if ignore_case:
        ser = ser.str.lower()
    return ser

def explode_if_needed(ser: pd.Series, split_regex: str | None):
    if not split_regex:
        # treat the entire cell as a single term
        return ser
    # split and explode into one-term-per-row Series
    split = ser.str.split(split_regex, regex=True)
    exploded = split.explode().dropna().astype(str).str.strip()
    # drop empties caused by leading/trailing delimiters
    exploded = exploded[exploded != ""]
    return exploded

def count_exact(go_descriptions, terms_series: pd.Series) -> dict:
    # Pre-compute value counts for speed
    vc = terms_series.value_counts(dropna=False)
    # Make sure we return zeros for missing terms
    return {desc: int(vc.get(desc, 0)) for desc in go_descriptions}

def process(go_desc_file, bacterial_files, sep, header, col_idx,
            split_regex, ignore_case, collapse_spaces, out_path):
    go_descriptions = read_go_descriptions(go_desc_file, ignore_case, collapse_spaces)

    # Prepare output table
    out_df = pd.DataFrame({"GO Description": go_descriptions})

    for path in bacterial_files:
        col = load_go_column(path, sep, header, col_idx)
        col = normalize_series(col, ignore_case, collapse_spaces)
        tokens = explode_if_needed(col, split_regex)
        counts = count_exact(go_descriptions, tokens)
        out_df[os.path.basename(path)] = out_df["GO Description"].map(counts)

    out_df.to_csv(out_path, index=False)
    print(f"[OK] Saved exact-match counts to {out_path}")
    print(f"[OK] Files processed: {len(bacterial_files)} | Terms listed: {len(go_descriptions)}")
    if split_regex:
        print(f"[OK] Split regex used: {split_regex}")

def parse_args():
    ap = argparse.ArgumentParser(
        description="Count GO descriptions by exact match (no substring/regex)."
    )
    ap.add_argument("go_descriptions", help="Text file with one GO description per line.")
    ap.add_argument("bacterial_files", nargs="+", help="TSV/CSV files with a GO column to count from.")
    ap.add_argument("--sep", default="\t", help="Input delimiter (default: tab)")
    ap.add_argument("--header", type=int, default=-1,
                    help="Header row index (default: -1 means no header). Use 0 if the first row is a header.")
    ap.add_argument("--col", type=int, default=3,
                    help="Zero-based index of the GO column (default: 3 => 4th column).")
    ap.add_argument("--split-regex", default=None,
                    help="Regex to split multiple terms per cell (e.g., '\\s*[;|,]\\s*'). If omitted, no split.")
    ap.add_argument("--ignore-case", action="store_true",
                    help="Case-insensitive matching (default: off).")
    ap.add_argument("--no-collapse-spaces", action="store_true",
                    help="Do not collapse internal whitespace (default: collapse).")
    ap.add_argument("--out", default="go_counts_output.csv",
                    help="Output CSV path (default: go_counts_output.csv)")
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    try:
        process(
            go_desc_file=args.go_descriptions,
            bacterial_files=args.bacterial_files,
            sep=args.sep,
            header=args.header,
            col_idx=args.col,
            split_regex=args.split_regex,
            ignore_case=args.ignore_case,
            collapse_spaces=not args.no_collapse_spaces,
            out_path=args.out,
        )
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
