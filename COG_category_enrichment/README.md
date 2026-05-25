# COG functional category enrichment analysis

This subdirectory contains the scripts and protocol for the COG functional category enrichment analysis described in §2.10.2 of Bisenić et al. (2026), corresponding to Figure 2C.

## Software requirements

- Python ≥3.8
- Python packages: `pandas`, `numpy`, `scipy`, `statsmodels`
- Microsoft Excel, for input file preparation and assembly of the combined count table
- GraphPad Prism v10, for the final bar plot

Scripts and commands were tested on Ubuntu inside Windows (WSL2).

## Files in this directory

- `cog_analysis.py` — assigns annotated genes to COG functional categories and produces per-genome gene-count tables
- `cog_enrichment_fisher.py` — Fisher's exact test for over-representation of COG categories in a favored group of genomes

## Workflow

### Step 1. Prepare per-genome COG input files

For each genome of interest:

1. Open the eggNOG-mapper annotation Excel file for the genome.
2. Copy columns **G** (`COG_category`), **H** (`Description`), and **L** (`KEGG_ko`) — in this order — into a new Excel workbook.
3. Delete the first three header rows.
4. With the new sheet selected, open Find & Replace (Ctrl + H). Enter `ko:` in "Find what", leave "Replace with" empty, and click "Replace All".
5. Save the resulting three-column table as a tab-delimited plain text file named `cog_genes_input.txt`.

### Step 2. Assign genes to COG categories

In the same working directory as `cog_genes_input.txt`, run:

```
python cog_analysis.py
```

The script reads `cog_genes_input.txt`, classifies each annotated gene into one or more of the standard top-level COG functional categories (genes assigned multiple COG letters are counted once per category), and writes:

- `cog_enrichment_results.csv` — gene counts per COG category for the current genome.
- `grouped_genes_by_cog.csv` — the genes themselves grouped by category, for inspection only.

**Repeat steps 1 and 2 for every genome to be compared.** Because `cog_analysis.py` reads and writes hardcoded filenames, before running the script for the next genome, rename `cog_enrichment_results.csv` to a genome-specific name (for example, `cog_enrichment_NGB245.csv`) so that it is not overwritten on the next run. Replace `cog_genes_input.txt` with the next genome's input before each run.

### Step 3. Assemble the combined count table

Build a single Excel workbook named `cog_enrichment_input.xlsx` combining the per-genome count tables:

- The first column, headed `COG Category`, contains the COG category letters (A, B, C, …, Z) in the same order they appear in the per-genome `cog_enrichment_*.csv` files. This can be copied directly from the first column of any one of those files.
- Each subsequent column corresponds to one of the 17 analyzed genomes. The column header is the genome name, and the cell values are the gene counts copied from the `Gene Count` column of that genome's `cog_enrichment_*.csv` file.

The genome names used as column headers must exactly match the names passed to `cog_enrichment_fisher.py` in step 4.

### Step 4. Test for COG category enrichment

Fisher's exact test for over-representation is run twice — once with the *Pilosibacter* genomes as the favored set, and once with the saccharolytic butyrogen reference genomes as the favored set. The script pools gene counts across the favored genomes and compares against the pooled gene counts across the remaining genomes, applying a one-sided Fisher's exact test (over-representation) followed by Benjamini–Hochberg correction across all COG categories.

To identify COG categories enriched in *Pilosibacter*:

```
python cog_enrichment_fisher.py cog_enrichment_input.xlsx \
  NGB245 P_fragilis P_sp900066055 P_sp902496665 P_sp934370915 P_sp963622995 P_sp963646925
```

Rename the resulting `cog_enrichment_fisher.csv` to `cog_enrichment_pilosibacter_favored.csv` before the next run, so that it is not overwritten.

To identify COG categories enriched in saccharolytic butyrogens:

```
python cog_enrichment_fisher.py cog_enrichment_input.xlsx \
  A_hadrus A_hallii A_rectalis C_catus C_eutactus F_prausnitzii R_hominis R_intestinalis R_inulinivorans S_variabile
```

Rename the resulting output to `cog_enrichment_saccharolytic_favored.csv`.

Each output CSV lists, per COG category: contingency table entries, group totals, group proportions, the odds ratio, the raw p-value, and the BH q-value. Categories are sorted by ascending q-value (most significant first).

Genome accessions used for the two groups are listed in Supplementary Table 2 of Bisenić et al. (2026).

### Step 5. Bar plot (Figure 2C)

COG categories reaching q ≤ 0.05 in either Fisher run were considered significantly enriched. For each significant category, per-genome gene counts (taken from `cog_enrichment_input.xlsx`) were normalized to gene counts per 1 Mb using the genome sizes reported in the Prokka assembly statistics for each genome. The resulting normalized values were imported into GraphPad Prism v10 to generate the bar plot shown in Figure 2C, with bars representing group means and individual data points representing per-genome values.
