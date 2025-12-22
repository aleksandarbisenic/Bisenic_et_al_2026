#!/usr/bin/env python3
# cog_enrichment_fisher.py
# Usage:
#   python cog_enrichment_fisher.py input.xlsx Favored1 [Favored2 ...]
#
# Input table format:
#   - First column: "COG Category"
#   - Remaining columns: genome names
#   - Cell values: INTEGER COUNTS (# of genes annotated to that COG in that genome)
#
# Assumptions:
#   * Each gene is assigned to at most one COG category (top-level categories).
#   * The table includes all COG-annotated genes (so totals = sum across COGs per genome).
#
# Test per COG:
#   Build a 2x2 on pooled gene counts across favored vs other genomes:
#       a = genes in favored group with this COG
#       b = total COG-annotated genes in favored group - a
#       c = genes in other group with this COG
#       d = total COG-annotated genes in other group - c
#   Fisher's exact test (alternative='greater') + BH-FDR over all COGs.
#
# Output:
#   cog_enrichment_fisher.csv with columns:
#     COG_Category, a_in_set, b_in_set_not, c_out_has, d_out_not,
#     total_in_set, total_out, prop_in_set, prop_out, odds_ratio, p_value, q_BH, n_favored, n_others
#
# Notes:
#   - If any genome column contains non-integers, the script will error.
#   - If totals differ from the sum of COGs (e.g., due to missing categories), provide a complete table.

import sys
import pandas as pd
import numpy as np

try:
    from scipy.stats import fisher_exact
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

try:
    from statsmodels.stats.multitest import multipletests
    _HAS_SM = True
except Exception:
    _HAS_SM = False

def load_any(path: str):
    pl = path.lower()
    if pl.endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(path, sheet_name=0)
    elif pl.endswith(".csv"):
        return pd.read_csv(path)
    else:
        raise ValueError("Unsupported file type. Use .xlsx/.xls/.csv")

def main():
    if len(sys.argv) < 3:
        print("Usage: python cog_enrichment_fisher.py input.xlsx Favored1 [Favored2 ...]")
        sys.exit(1)
    input_path = sys.argv[1]
    favored = sys.argv[2:]

    df = load_any(input_path)
    df.columns = [str(c) for c in df.columns]
    if "COG Category" not in df.columns or df.shape[1] < 3:
        raise SystemExit("Expected first column 'COG Category' followed by genome columns.")
    genomes = [c for c in df.columns if c != "COG Category"]

    # Validate favored present
    missing = [g for g in favored if g not in genomes]
    if missing:
        raise SystemExit(f"Favored genomes not found: {missing}. Present: {genomes}")

    # Ensure integer counts
    for g in genomes:
        if not np.issubdtype(df[g].dropna().dtype, np.integer):
            # Try casting safely
            vals = df[g].astype(float)
            if not np.all(np.isclose(vals, np.round(vals))):
                raise SystemExit(f"Column '{g}' contains non-integer values. Provide integer gene counts per COG.")
            df[g] = vals.astype(int)

    # Totals per group (sum of all COG counts)
    favored_set = set(favored)
    in_cols = [g for g in genomes if g in favored_set]
    out_cols = [g for g in genomes if g not in favored_set]

    total_in = int(df[in_cols].sum().sum())
    total_out = int(df[out_cols].sum().sum())

    records = []
    for _, row in df.iterrows():
        cog = row["COG Category"]
        a = int(row[in_cols].sum()) if in_cols else 0
        c = int(row[out_cols].sum()) if out_cols else 0
        b = total_in - a
        d = total_out - c
        table = np.array([[a, b], [c, d]], dtype=int)

        if _HAS_SCIPY:
            odds, p = fisher_exact(table, alternative="greater")
        else:
            odds = (a*d) / max(1, b*c)
            p = np.nan

        prop_in = a / total_in if total_in > 0 else np.nan
        prop_out = c / total_out if total_out > 0 else np.nan

        records.append({
            "COG_Category": cog,
            "a_in_set": a,
            "b_in_set_not": b,
            "c_out_has": c,
            "d_out_not": d,
            "total_in_set": total_in,
            "total_out": total_out,
            "prop_in_set": prop_in,
            "prop_out": prop_out,
            "odds_ratio": odds,
            "p_value": p,
            "n_favored": len(in_cols),
            "n_others": len(out_cols)
        })

    out = pd.DataFrame(records)

    # BH-FDR
    if out["p_value"].notna().any():
        if _HAS_SM:
            out["q_BH"] = multipletests(out["p_value"].to_numpy(), method="fdr_bh")[1]
        else:
            p = out["p_value"].to_numpy()
            order = np.argsort(p)
            ranks = np.empty_like(order); ranks[order] = np.arange(1, len(p)+1)
            q = p * len(p) / ranks
            q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
            q_out = np.empty_like(q_sorted); q_out[order] = q_sorted
            out["q_BH"] = np.clip(q_out, 0, 1)
    else:
        out["q_BH"] = np.nan

    out = out.sort_values(["q_BH","p_value","odds_ratio"], ascending=[True, True, False])
    out.to_csv("cog_enrichment_fisher.csv", index=False)

    print(f"[OK] Processed {len(out)} COG categories across {len(genomes)} genomes.")
    print(f"[OK] Favored group total COG-annotated genes: {total_in}; others: {total_out}.")
    if _HAS_SCIPY:
        print("[OK] Fisher exact (one-sided) run; BH-FDR applied. Results -> cog_enrichment_fisher.csv")
    else:
        print("[WARN] scipy not found; p-values are NaN. Install scipy to compute exact p-values.")

if __name__ == "__main__":
    main()
