import pandas as pd

# Load your binary and percentage CSV files
binary_df = pd.read_csv('kegg_module_completeness_binary.csv', index_col=0)
percent_df = pd.read_csv('kegg_module_completeness_percentage.csv', index_col=0)

# Identify columns to keep: those with at least one "1" and at least one "0" in the binary file
columns_to_keep = binary_df.columns[(binary_df == 1).any(axis=0) & (binary_df == 0).any(axis=0)]

# Filter both DataFrames to keep only the identified columns
binary_df_filtered = binary_df[columns_to_keep]
percent_df_filtered = percent_df[columns_to_keep]

# Save the filtered DataFrames to new CSV files
binary_df_filtered.to_csv('differentially_present_binary.csv')
percent_df_filtered.to_csv('differentially_present_percentage.csv')

print("Filtered data has been saved to 'kegg_module_completeness_binary_filtered.csv' and 'kegg_module_completeness_percentage_filtered.csv'")
