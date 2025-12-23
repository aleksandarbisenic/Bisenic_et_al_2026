#!/usr/bin/env python3
import sys
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

if len(sys.argv) != 3:
    sys.exit("Usage: python plot_upset.py <matrix.tsv> <GenomeColumn>")

matrix_tsv, target = sys.argv[1:3]
df = pd.read_csv(matrix_tsv, sep="\t").astype(bool)
if target not in df.columns:
    sys.exit(f"Column '{target}' not found in {list(df.columns)}")

# 1) full Series of all intersections
series_full = from_indicators(df.columns.tolist(), df)

# 2) plot with only the top-10 intersections drawn
plt.figure(figsize=(10, 6))
UpSet(
    series_full,
    sort_by="cardinality",
    max_subset_rank=10,   # ← new in v0.9.0 to show only N largest bars
    show_counts="%d"
).plot()

plt.suptitle("Top-10 gene-family intersections across seven Pilosibacter genomes",
             fontsize=14)
plt.tight_layout()
plt.savefig("upset_pilosibacter_top10.svg")
plt.savefig("upset_pilosibacter_top10.png", dpi=300)
print("Saved upset_pilosibacter_top10.svg / .png")
