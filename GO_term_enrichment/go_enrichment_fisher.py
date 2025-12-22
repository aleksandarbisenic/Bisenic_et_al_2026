#!/usr/bin/env python3
# go_enrichment_fisher.py
# Usage: python go_enrichment_fisher.py input_file.xlsx NGB244 NGB245 [more favored genomes...]
# Input can be:
#   (A) Wide paired columns (like: NGB244, "GO term", NGB241, "GO term.1", ...)
#   (B) Long format with columns: genome, gene_id, go_id(or "GO term")
# Output: enriched_GO.csv with Fisher one-sided (greater) enrichment for GO terms in favored genomes vs others.

import argparse
import sys
import pandas as pd
import numpy as np

# Try to import stats packages; fall back gracefully if missing
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

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR (fallback if statsmodels is missing)."""
    p = np.asarray(pvals, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranks = np.empty(m, dtype=int); ranks[order] = np.arange(1, m+1)
    q = p * m / ranks
    q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
    out = np.empty_like(q_sorted); out[order] = q_sorted
    return np.clip(out, 0, 1)

def load_any(path: str, sheet: str | None) -> pd.DataFrame:
    pl = path.lower()
    if pl.endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(path, sheet_name=sheet or 0)
    elif pl.endswith(".csv") or pl.endswith(".tsv"):
        sep = "," if pl.endswith(".csv") else "\t"
        return pd.read_csv(path, sep=sep)
    else:
        raise ValueError("Unsupported file type. Use .xlsx/.xls/.csv/.tsv")

def detect_format(df: pd.DataFrame) -> str:
    cols = df.columns.str.lower().tolist()
    if {"genome","gene_id"}.issubset(set(cols)) and ("go_id" in cols or "go term" in cols):
        return "long"
    # Wide paired if we see multiple "go term" columns interleaved
    go_like = [c for c in df.columns if str(c).lower().startswith("go term")]
    if go_like and len(go_like) >= 1:
        return "wide"
    # Heuristic: two-column repeated pattern
    if df.shape[1] >= 2:
        return "wide"
    return "unknown"

def to_long_from_wide(df: pd.DataFrame) -> pd.DataFrame:
    """Convert wide paired columns into long format with columns: genome, gene_id, go_term."""
    cols = list(df.columns)
    long_rows = []
    i = 0
    # Strategy: pair each non-GO column with the immediate next GO-like column
    while i < len(cols) - 1:
        gene_col = cols[i]
        go_col = cols[i+1]
        if str(go_col).lower().startswith("go term"):
            genome = str(gene_col)
            sub = df[[gene_col, go_col]].dropna(subset=[gene_col, go_col])
            sub.columns = ["gene_id", "go_term"]
            sub["genome"] = genome
            long_rows.append(sub[["genome","gene_id","go_term"]])
            i += 2
        else:
            # skip this column if it doesn't have a GO partner
            i += 1
    if not long_rows:
        raise ValueError("Could not parse paired GO columns. Ensure the file alternates <genome>, 'GO term' columns.")
    out = pd.concat(long_rows, ignore_index=True)
    # Deduplicate identical triplets
    out = out.drop_duplicates(subset=["genome","gene_id","go_term"])
    return out

def to_long(df: pd.DataFrame) -> pd.DataFrame:
    fmt = detect_format(df)
    if fmt == "long":
        cols = {c.lower(): c for c in df.columns}
        genome_col = cols["genome"]
        gene_col = cols["gene_id"]
        go_col = cols["go_id"] if "go_id" in cols else cols.get("go term")
        if go_col is None:
            raise ValueError("Long format requires a 'go_id' or 'GO term' column.")
        out = df[[genome_col, gene_col, go_col]].copy()
        out.columns = ["genome","gene_id","go_term"]
        # Drop NA and duplicates
        out = out.dropna(subset=["genome","gene_id","go_term"]).drop_duplicates()
        return out
    elif fmt == "wide":
        return to_long_from_wide(df)
    else:
        raise ValueError("Unrecognized input format. Provide either long format or wide paired columns.")

def build_counts(long_df: pd.DataFrame, favored: list[str]) -> pd.DataFrame:
    """Construct 2x2 counts for each GO term; items are unique (genome,gene_id) pairs."""
    # Clean strings
    long_df["genome"] = long_df["genome"].astype(str)
    long_df["gene_id"] = long_df["gene_id"].astype(str)
    long_df["go_term"] = long_df["go_term"].astype(str)

    # Universe: all unique genes that have any GO term
    # Define "gene key" as (genome, gene_id)
    long_df = long_df.drop_duplicates(subset=["genome","gene_id","go_term"])
    all_genes = long_df[["genome","gene_id"]].drop_duplicates()
    all_genes["in_set"] = all_genes["genome"].isin(favored)

    # For each GO, which genes have it?
    # Build mapping gene_key -> set of GO terms
    # But we can compute counts via groupby for speed
    has_go = long_df.assign(in_set=long_df["genome"].isin(favored))
    # Count a: in_set & has_go; c: out_set & has_go
    go_counts = has_go.groupby(["go_term","in_set"], as_index=False)[["gene_id"]].nunique()
    # Pivot to columns: in_set True/False
    pivot = go_counts.pivot(index="go_term", columns="in_set", values="gene_id").fillna(0).astype(int)
    pivot = pivot.rename(columns={True:"a_in_set_has", False:"c_out_has"})

    # Compute totals for set sizes
    n_set = int(all_genes["in_set"].sum())       # genes in favored genomes
    n_out = int((~all_genes["in_set"]).sum())    # genes in other genomes

    # Derive b and d
    pivot["b_in_set_not"] = n_set - pivot["a_in_set_has"].astype(int)
    pivot["d_out_not"]   = n_out - pivot["c_out_has"].astype(int)

    pivot = pivot[["a_in_set_has","b_in_set_not","c_out_has","d_out_not"]].reset_index()
    return pivot, n_set, n_out

def run_fisher(df_counts: pd.DataFrame) -> pd.DataFrame:
    """Run Fisher one-sided (greater) enrichment per GO term."""
    pvals = []
    ors = []
    for _, row in df_counts.iterrows():
        a,b,c,d = int(row["a_in_set_has"]), int(row["b_in_set_not"]), int(row["c_out_has"]), int(row["d_out_not"])
        table = np.array([[a,b],[c,d]], dtype=int)
        if _HAS_SCIPY:
            # scipy returns oddsratio and p-value
            # Use alternative='greater' to test enrichment in favored set
            odds, p = fisher_exact(table, alternative="greater")
        else:
            # Simple approximation if scipy unavailable
            odds = (a*d) / max(1, b*c)
            # No exact p without scipy; use placeholder
            # (User should install scipy for exact p-values)
            p = np.nan
        pvals.append(p); ors.append(odds)
    out = df_counts.copy()
    out["odds_ratio"] = ors
    out["p_value"] = pvals
    # FDR
    p_arr = np.array(pvals, dtype=float)
    if np.isnan(p_arr).all():
        out["q_BH"] = np.nan
    else:
        if _HAS_SM:
            out["q_BH"] = multipletests(p_arr, method="fdr_bh")[1]
        else:
            out["q_BH"] = bh_fdr(p_arr)
    return out

def attach_support(long_df: pd.DataFrame, favored: list[str], results: pd.DataFrame, max_list: int = 20) -> pd.DataFrame:
    """Add support counts and (optionally truncated) gene lists for favored genomes."""
    # Build mapping GO -> list of genes in favored set that have it
    mask = long_df["genome"].isin(favored)
    fav = long_df[mask][["genome","gene_id","go_term"]].drop_duplicates()
    grp = fav.groupby("go_term")
    support_counts = grp.size().rename("genes_in_set_with_GO").reset_index()
    # Build short lists "genome:gene" up to max_list
    lists = grp.apply(lambda g: ";".join((g["genome"] + ":" + g["gene_id"]).head(max_list))).rename("example_genes").reset_index()
    out = results.merge(support_counts, on="go_term", how="left").merge(lists, on="go_term", how="left")
    out["genes_in_set_with_GO"] = out["genes_in_set_with_GO"].fillna(0).astype(int)
    out["example_genes"] = out["example_genes"].fillna("")
    return out

def main():
    ap = argparse.ArgumentParser(description="GO over-representation (Fisher) for favored genomes vs others.")
    ap.add_argument("input", help="Input Excel/CSV (wide paired columns or long format)")
    ap.add_argument("favored", nargs="+", help="Names of favored genomes (columns or 'genome' values)")
    ap.add_argument("--sheet", default=None, help="Excel sheet name (if Excel)")
    ap.add_argument("--out", default="enriched_GO.csv", help="Output CSV (default enriched_GO.csv)")
    args = ap.parse_args()

    # Load
    df = load_any(args.input, args.sheet)

    # Convert to long
    long_df = to_long(df)

    # Validate favored names present
    present = set(long_df["genome"].unique().tolist())
    missing = [g for g in args.favored if g not in present]
    if missing:
        sys.exit(f"Favored genome(s) not found: {missing}. Present genomes: {sorted(present)}")

    # Build 2x2 counts
    counts_df, n_set, n_out = build_counts(long_df, args.favored)

    # Run Fisher
    res = run_fisher(counts_df)

    # Attach support info
    res = attach_support(long_df, args.favored, res)

    # Sort by q then p then odds
    res = res.sort_values(["q_BH","p_value","odds_ratio"], ascending=[True, True, False])

    # Write
    res.to_csv(args.out, index=False)

    # Console summary
    n_terms = res.shape[0]
    n_sig = int((res["q_BH"] < 0.05).sum()) if "q_BH" in res.columns else 0
    print(f"[OK] Processed {n_terms} GO terms across {len(present)} genomes.")
    print(f"[OK] Favored set size (genes): {n_set}; others: {n_out}.")
    if _HAS_SCIPY:
        print(f"[OK] Fisher exact (one-sided) run; BH-FDR applied. Significant (q<0.05): {n_sig}.")
    else:
        print("[WARN] scipy not found; p-values are NaN. Install scipy for exact Fisher p-values.")
    print(f"[OK] Results written to {args.out}")

if __name__ == "__main__":
    main()
