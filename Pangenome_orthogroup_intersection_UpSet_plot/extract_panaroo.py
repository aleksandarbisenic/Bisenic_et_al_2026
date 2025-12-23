#!/usr/bin/env python3
"""
Create:
  (1) Boolean presence/absence matrix
  (2) list of gene families unique to a chosen genome (with locus IDs)
  (3) list of core gene families present in ALL genomes (with locus IDs)

Usage:
  python panaroo_parse.py <gene_presence_absence.csv> <GenomeColumn>

Example:
  python panaroo_parse.py gene_presence_absence.csv 245_genome
"""

import sys
import pandas as pd

if len(sys.argv) != 3:
    sys.exit("Usage: python panaroo_parse.py <gene_presence_absence.csv> <GenomeColumn>")

CSV_IN  = sys.argv[1]
TARGET  = sys.argv[2]           # e.g. "245_genome"

MATRIX_OUT = "upset_matrix.tsv"
UNIQUE_OUT = f"{TARGET}_unique_genes.tsv"
CORE_OUT   = "core_genes.tsv"

# 1. Read the CSV; keep everything as string so blanks stay NaN
df = pd.read_csv(CSV_IN, dtype=str)

# 2. Identify genome columns: everything after the first 3 metadata cols
meta_cols   = list(df.columns[:3])                # ["Gene","Non-unique Gene name","Annotation"]
genome_cols = [c for c in df.columns if c not in meta_cols]

if TARGET not in genome_cols:
    sys.exit(f"ERROR: column '{TARGET}' not found.\nAvailable genome columns:\n" +
             ", ".join(genome_cols))

# 3. Build Boolean presence/absence matrix (1 = present, 0 = absent)
bool_df = df[genome_cols] \
            .fillna("") \
            .applymap(lambda x: 1 if x.strip() else 0)

bool_df.to_csv(MATRIX_OUT, sep="\t", index=False)
print(f"[OK] Boolean matrix written to {MATRIX_OUT}")

# 4. Extract genes unique to TARGET
others = [c for c in genome_cols if c != TARGET]
mask_unique = (bool_df[TARGET] == 1) & (bool_df[others].sum(axis=1) == 0)

unique_df = df.loc[mask_unique, meta_cols + [TARGET]].copy()
unique_df = unique_df.rename(columns={TARGET: f"Locus_in_{TARGET}"})
unique_df.to_csv(UNIQUE_OUT, sep="\t", index=False)
print(f"[OK] {len(unique_df)} families unique to {TARGET} saved to {UNIQUE_OUT}")

# 5. Extract core genes (present in *all* genomes)
mask_core = (bool_df.sum(axis=1) == len(genome_cols))

# keep all locus columns for core genes
core_df = df.loc[mask_core, meta_cols + genome_cols].copy()
# rename each genome column to include "Locus_in_"
core_df = core_df.rename(columns={g: f"Locus_in_{g}" for g in genome_cols})
core_df.to_csv(CORE_OUT, sep="\t", index=False)
print(f"[OK] {len(core_df)} core families saved to {CORE_OUT}")
