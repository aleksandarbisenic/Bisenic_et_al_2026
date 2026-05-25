#!/usr/bin/env python3
"""
kegg_module_fisher_transposed.py

✅ This version is for the *TRANSPOSED* format:

- Rows = genomes
- Columns = KEGG modules
- Cells = 0/1 (complete/incomplete)

Example:

Genome              ModuleA   ModuleB   ModuleC
A_hadrus            1         0         1
A_hallii            1         1         1
...
P_sp963646925       0         1         0

It runs Fisher's exact test per module (2×2: complete vs incomplete × group A vs group B),
then Benjamini–Hochberg FDR correction across modules.

USAGE (Excel/CSV):
python kegg_module_fisher_transposed.py --in kegg_binary.xlsx --out fisher_results.csv \
  --groupA A_hadrus A_hallii A_rectalis C_catus C_eutactus F_prausnitzii R_hominis R_intestinalis R_inulinivorans S_variabile \
  --groupB NGB245 P_fragilis P_sp900066055 P_sp902496665 P_sp934370915 P_sp963622995 P_sp963646925
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

try:
    from statsmodels.stats.multitest import multipletests
except ImportError:
    print("ERROR: statsmodels is required. Install with: pip install statsmodels", file=sys.stderr)
    sys.exit(1)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR q-values."""
    q = np.full_like(pvals, np.nan, dtype=float)
    mask = np.isfinite(pvals)
    if mask.sum() == 0:
        return q
    _, qvals, _, _ = multipletests(pvals[mask], alpha=0.05, method="fdr_bh")
    q[mask] = qvals
    return q


def read_input_table(infile: str, sep: str | None, sheet: str | None) -> pd.DataFrame:
    """Read CSV/TSV or Excel based on extension."""
    path = Path(infile)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {infile}")

    ext = path.suffix.lower()
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(infile, sheet_name=sheet)
    return pd.read_csv(infile, sep=sep, engine="python")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="infile", required=True, help="Input CSV/TSV/Excel file.")
    ap.add_argument("--out", dest="outfile", required=True, help="Output CSV file.")
    ap.add_argument("--sep", default=None, help="Delimiter for CSV/TSV (optional).")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (optional).")

    ap.add_argument(
        "--genome_col",
        default=None,
        help="Column containing genome names. If omitted, uses the FIRST column.",
    )

    ap.add_argument("--groupA", nargs="*", default=None, help="Genome names for Group A.")
    ap.add_argument("--groupB", nargs="*", default=None, help="Genome names for Group B.")

    ap.add_argument(
        "--alternative",
        choices=["two-sided", "greater", "less"],
        default="two-sided",
        help="Fisher exact test alternative.",
    )

    ap.add_argument(
        "--min_total_complete",
        type=int,
        default=0,
        help="Only test modules with at least this many total 1s across all genomes.",
    )
    ap.add_argument(
        "--max_total_complete",
        type=int,
        default=None,
        help="Only test modules with at most this many total 1s across all genomes.",
    )

    return ap.parse_args()


def main():
    args = parse_args()

    # ---- Read table ----
    try:
        df = read_input_table(args.infile, args.sep, args.sheet)
    except Exception as e:
        print(f"ERROR reading file: {e}", file=sys.stderr)
        sys.exit(1)

    if df.shape[1] < 3:
        print("ERROR: Need >= 1 genome column + >= 2 module columns.", file=sys.stderr)
        sys.exit(1)

    genome_col = args.genome_col if args.genome_col else df.columns[0]
    if genome_col not in df.columns:
        print(f"ERROR: genome_col '{genome_col}' not found. Columns: {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    module_cols = [c for c in df.columns if c != genome_col]
    if len(module_cols) < 1:
        print("ERROR: No module columns detected.", file=sys.stderr)
        sys.exit(1)

    # ---- Set genome names as index ----
    work = df.copy()
    work[genome_col] = work[genome_col].astype(str)
    work = work.set_index(genome_col)

    # ---- Coerce module values to 0/1 ----
    for c in module_cols:
        work[c] = pd.to_numeric(work[c], errors="coerce").fillna(0)
        work[c] = (work[c] > 0).astype(int)

    # ---- Validate groups ----
    if not args.groupA or not args.groupB:
        print("ERROR: You must provide BOTH --groupA and --groupB genome name lists.", file=sys.stderr)
        sys.exit(1)

    missingA = [g for g in args.groupA if g not in work.index]
    missingB = [g for g in args.groupB if g not in work.index]
    if missingA or missingB:
        print("ERROR: Some genomes in groups are not found in the genome rows.", file=sys.stderr)
        if missingA:
            print(f"  Missing Group A genomes: {missingA}", file=sys.stderr)
        if missingB:
            print(f"  Missing Group B genomes: {missingB}", file=sys.stderr)
        sys.exit(1)

    groupA = work.loc[args.groupA, module_cols]
    groupB = work.loc[args.groupB, module_cols]
    nA = groupA.shape[0]
    nB = groupB.shape[0]

    # ---- Fisher per module ----
    rows = []
    for module in module_cols:
        a = int(groupA[module].sum())   # complete in A
        b = int(groupB[module].sum())   # complete in B

        total_complete = a + b
        if total_complete < args.min_total_complete:
            continue
        if args.max_total_complete is not None and total_complete > args.max_total_complete:
            continue

        table = [[a, nA - a], [b, nB - b]]
        or_val, p = fisher_exact(table, alternative=args.alternative)

        prevA = a / nA
        prevB = b / nB
        delta = prevA - prevB

        rows.append(
            {
                "Module": module,
                "GroupA_complete": a,
                "GroupA_total": nA,
                "GroupA_prevalence": prevA,
                "GroupB_complete": b,
                "GroupB_total": nB,
                "GroupB_prevalence": prevB,
                "Delta_prevalence_A_minus_B": delta,
                "Odds_ratio": or_val,
                "P_value": p,
            }
        )

    if not rows:
        print("No modules passed filters.", file=sys.stderr)
        sys.exit(1)

    res = pd.DataFrame(rows)
    res["Q_value_BH_FDR"] = bh_fdr(res["P_value"].values)

    res = res.sort_values(["Q_value_BH_FDR", "P_value", "Module"], ascending=[True, True, True])
    res.to_csv(args.outfile, index=False)

    sig = res[res["Q_value_BH_FDR"] <= 0.05]
    print("Done.")
    print(f"Tested modules: {len(res)}")
    print(f"Significant (q<=0.05): {len(sig)}")

    if len(sig) > 0:
        print("\nTop significant modules:")
        show = sig.head(15)[
            [
                "Module",
                "GroupA_complete",
                "GroupB_complete",
                "Delta_prevalence_A_minus_B",
                "Odds_ratio",
                "P_value",
                "Q_value_BH_FDR",
            ]
        ]
        with pd.option_context("display.max_colwidth", 80, "display.width", 140):
            print(show.to_string(index=False))


if __name__ == "__main__":
    main()
