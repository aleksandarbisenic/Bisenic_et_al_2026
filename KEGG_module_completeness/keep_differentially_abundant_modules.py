import pandas as pd

# Load your binary and percentage CSV files
binary_df = pd.read_csv('kegg_module_completeness_binary.csv', index_col=0)
percent_df = pd.read_csv('kegg_module_completeness_percentage.csv', index_col=0)

# Sanity check: both files should have the same rows (genomes) and columns (modules)
if binary_df.shape != percent_df.shape:
    raise ValueError(
        "Binary and percentage tables have different shapes: "
        f"{binary_df.shape} vs {percent_df.shape}. "
        "Make sure they contain the same genomes (rows) and modules (columns)."
    )

n_genomes = binary_df.shape[0]

# Base rule (previous behavior):
# Keep modules that are not invariant (i.e., at least one 1 and at least one 0).
has_one = (binary_df == 1).any(axis=0)
has_zero = (binary_df == 0).any(axis=0)
keep_mask = has_one & has_zero

# New rule:
# If there are >10 genomes, also remove modules that differ in only one genome:
# - exactly one 1 (all others 0), OR
# - exactly one 0 (all others 1).
if n_genomes > 10:
    ones_count = (binary_df == 1).sum(axis=0)
    zeros_count = (binary_df == 0).sum(axis=0)
    keep_mask = keep_mask & ~((ones_count == 1) | (zeros_count == 1))

columns_to_keep = binary_df.columns[keep_mask]

# Filter both DataFrames to keep only the identified columns
binary_df_filtered = binary_df.loc[:, columns_to_keep]
percent_df_filtered = percent_df.loc[:, columns_to_keep]

# Save the filtered DataFrames to new CSV files
binary_out = 'differentially_present_binary.csv'
percent_out = 'differentially_present_percentage.csv'
binary_df_filtered.to_csv(binary_out)
percent_df_filtered.to_csv(percent_out)

dropped = binary_df.shape[1] - len(columns_to_keep)
print(
    f"Genomes: {n_genomes}\n"
    f"Kept modules: {len(columns_to_keep)} / {binary_df.shape[1]} (dropped {dropped})\n"
    "Saved to:\n"
    f"  - {binary_out}\n"
    f"  - {percent_out}"
)
