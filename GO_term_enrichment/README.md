# GO term enrichment analysis

This subdirectory contains the scripts and protocol for the GO term enrichment analysis described in §2.10.3 of Bisenić et al. (2026), corresponding to Figures 2D and 2E.

## Software requirements

- Python ≥3.8
- Python packages: `pandas`, `numpy`, `scipy`, `statsmodels`, `openpyxl`
- Microsoft Excel with VBA enabled, for filtering the PANNZER2 GO files to Biological Process terms
- GraphPad Prism v10, for the final heatmaps

Scripts and commands were tested on Ubuntu inside Windows (WSL2).

## Files in this directory

- `pannzer2_GO_file_parser_script.txt` — Excel VBA macro for filtering the PANNZER2 GO output to Biological Process (BP) ontology rows
- `go_enrichment_fisher.py` — Fisher's exact test for over-representation of GO terms in a favored group of genomes
- `differential_go_terms.py` — counts the number of unique genes annotated with each GO description in each genome
- `go_log1p_conversion.py` — applies natural log1p transformation to gene count tables for heatmap visualization
- `go_term_to_gene.py` — extracts gene-level evidence for each significant GO term from the PANNZER2 annotation files (exact GO term match)
- `level2_gene_search.py` — substring search within gene description rows only, case-insensitive
- `level3_gene_search.py` — substring search across all annotation columns, case-sensitive by default

## Workflow

### Step 1. Generate PANNZER2 annotations

For each genome of interest, submit its Prokka-generated `.faa` protein file to PANNZER2 (http://ekhidna2.biocenter.helsinki.fi/sanspanz/, accessed 21 January 2025) and download both output files:

- The **GO file**, used in steps 2–7 below for the Fisher enrichment analysis and heatmap visualization. It contains one row per gene–GO association, with the gene ID in column 1, the ontology (`BP`, `CC`, or `MF`) in column 2, and the GO term description in column 4. Save this file in Excel format, named after the genome (for example, `NGB245_GO.xlsx`).
- The **combined annotation file**, used in step 9 for the targeted gene searches. It contains multiple rows per gene with different row types in column B (`original_DE`, `DE`, `GN`, `BP_ARGOT`, `CC_ARGOT`, `MF_ARGOT`, etc.), with the corresponding text in column F. Save this file in Excel format in a separate folder dedicated to the targeted searches (for example, `pannzer2_combined/NGB245.xlsx`).

Both file types must be preserved — the macro and Fisher pipeline only need the GO file, but the targeted searches in step 9 require the combined annotation file because the BP-filtering macro in step 2 destroys the row types they depend on.

### Step 2. Filter the GO file to Biological Process (BP) ontology

For the enrichment analysis, only Biological Process annotations are retained from each genome's GO file.

For each genome's GO file:

1. Make a copy of the GO file (e.g., `NGB245_GO_BP.xlsx`) so that the original is preserved as a fallback.
2. Open the copy in Excel.
3. Open the Visual Basic editor (Alt + F11). From the menu, select Insert > Module. Paste the entire contents of `pannzer2_GO_file_parser_script.txt` into the module window, then close the editor.
4. Run the macro by pressing Alt + F8, selecting `FilterBP`, and clicking Run. The macro deletes all rows in which column B is not `BP`, leaving only Biological Process annotations.
5. Save the filtered workbook.

Repeat for every genome.

### Step 3. Assemble Fisher_GO_input.xlsx

Build a single Excel workbook named `Fisher_GO_input.xlsx` combining the BP-filtered annotations of all 17 analyzed genomes:

- For each genome, copy columns 1 and 4 of its BP-filtered workbook (gene ID and GO term description, respectively).
- Paste them as a pair of columns into the combined workbook.
- Rename the header of each gene-ID column to the genome name (e.g., `NGB245`).
- Name the first GO term column `GO term`. Subsequent ones are renamed automatically by Excel to `GO term.1`, `GO term.2`, and so on as you paste them in.

Each genome therefore occupies two adjacent columns in the combined workbook. The genome names used as headers must exactly match the names passed to `go_enrichment_fisher.py` in step 4.

### Step 4. Test for enriched GO terms

Fisher's exact test for over-representation is run twice — once with the *Pilosibacter* genomes as the favored set, and once with the saccharolytic butyrogen reference genomes as the favored set. The script pools unique-gene evidence across the favored genomes and compares against the pooled evidence across the remaining genomes, applying a one-sided Fisher's exact test followed by Benjamini–Hochberg correction across all GO terms.

To identify GO terms enriched in *Pilosibacter*:

```
python go_enrichment_fisher.py Fisher_GO_input.xlsx \
  NGB245 P_fragilis P_sp900066055 P_sp902496665 P_sp934370915 P_sp963622995 P_sp963646925
```

Rename the resulting `enriched_GO.csv` to `enriched_GO_pilosibacter_favored.csv` before the next run, so that it is not overwritten.

To identify GO terms enriched in saccharolytic butyrogens:

```
python go_enrichment_fisher.py Fisher_GO_input.xlsx \
  A_hadrus A_hallii A_rectalis C_catus C_eutactus F_prausnitzii R_hominis R_intestinalis R_inulinivorans S_variabile
```

Rename the resulting output to `enriched_GO_saccharolytic_favored.csv`.

Each output CSV lists, per GO term: contingency table entries, the odds ratio, the raw p-value, the BH q-value, and supporting gene examples from the favored set. GO terms are sorted by ascending q-value (most significant first).

Genome accessions used for the two groups are listed in Supplementary Table 2 of Bisenić et al. (2026).

### Step 5. Extract significant GO term descriptions

For each Fisher run, collect the GO term descriptions reaching q ≤ 0.05 into a plain text file, one description per line:

- `go_descriptions_pilosibacter.txt`
- `go_descriptions_saccharolytic.txt`

These lists are the inputs to the subsequent counting and gene-extraction steps.

### Step 6. Count genes per GO term per genome

For each direction of enrichment, run `differential_go_terms.py` against `Fisher_GO_input.xlsx`, using the corresponding significant-terms list:

```
python differential_go_terms.py go_descriptions_pilosibacter.txt Fisher_GO_input.xlsx
```

Rename the output `go_counts_output.csv` to `go_counts_pilosibacter.csv` before running the saccharolytic direction:

```
python differential_go_terms.py go_descriptions_saccharolytic.txt Fisher_GO_input.xlsx
```

Rename the resulting output to `go_counts_saccharolytic.csv`.

Each output is a count matrix with one row per significant GO term and one column per genome, where each cell is the number of unique genes in that genome annotated with that GO term.

### Step 7. Log-transform counts for heatmap visualization

Apply the natural log1p transformation (ln[count + 1]) to each count matrix:

```
python go_log1p_conversion.py go_counts_pilosibacter.csv
python go_log1p_conversion.py go_counts_saccharolytic.csv
```

The transformed tables are written to the `log_out/` subdirectory as `go_counts_pilosibacter_log1p.csv` and `go_counts_saccharolytic_log1p.csv`.

### Step 8. Heatmaps (Figures 2D and 2E)

The log1p-transformed count tables from step 7 were imported into GraphPad Prism v10 to generate the heatmaps shown in Figure 2D (GO terms enriched in *Pilosibacter*) and Figure 2E (GO terms enriched in saccharolytic butyrogens).

### Step 9. Targeted gene searches (*Pilosibacter*-enriched terms only)

The Fisher test in step 4 identifies GO terms over-represented in *Pilosibacter*, but a significant statistical result does not by itself confirm a coherent biological signal — it may reflect annotation idiosyncrasies (e.g., the same protein labeled with different terminology across genomes, or one gene assigned multiple closely related GO descriptions). To distinguish real biological signal from annotation artifacts, the genes underlying each Pilosibacter-enriched function were inspected directly across all 17 genomes using a tiered set of searches against the **combined PANNZER2 annotation files** preserved in step 1 (the files described in step 1 as `pannzer2_combined/<genome>.xlsx`).

This step was performed only for the *Pilosibacter*-enriched direction, as the genus was the functional focus of the manuscript.

The three search scripts implement progressively broader matching strategies, used together to investigate each enriched function thoroughly:

**Level 1 — exact GO term match (`go_term_to_gene.py`)**

Retrieves, from each genome's annotation file, the gene IDs and descriptions annotated with each GO term listed in a text file. Matching is on the exact GO term description.

```
python go_term_to_gene.py <pannzer2_combined_folder> go_descriptions_pilosibacter.txt
```

Use `go_descriptions_pilosibacter.txt` from step 5 as the input. Output: `go_hits_adjacent_columns.csv`, a single CSV with one block per GO term, each block containing paired columns (gene ID + description) for every genome. This identifies which genes account for the Fisher enrichment of each term in the Pilosibacter genomes and shows whether the term is absent or annotated differently in the saccharolytic comparators.

**Level 2 — substring match within description rows, case-insensitive (`level2_gene_search.py`)**

When a Level 1 result looks sparse or suspicious, this script searches across all genomes for gene-description substrings related to the function in question. Searches are restricted to `original_DE` rows (the original gene descriptions) and are case-insensitive.

```
python level2_gene_search.py <pannzer2_combined_folder> search_terms.txt
```

The `search_terms.txt` file is curated manually for each function being investigated — for example, when investigating biotin biosynthesis, the search terms might include the names of biotin pathway proteins (`Adenosylmethionine-8-amino-7-oxononanoate aminotransferase`, `ATP-dependent dethiobiotin synthetase`, and so on). Output: `term_hits_contains_originalDE.csv`, structured the same way as the Level 1 output.

**Level 3 — substring match across all annotation rows, case-sensitive (`level3_gene_search.py`)**

The broadest search. Looks for substring matches in any value in column F, regardless of row type, and is case-sensitive by default. Used when a function is suspected to be encoded in a genome but annotated under less standard terminology, where Level 2 returns nothing.

```
python level3_gene_search.py <pannzer2_combined_folder> search_terms.txt
```

Pass `--ignore-case` to switch to case-insensitive matching. Output: `term_hits_contains.csv`.

All three scripts read every `.xlsx` file in the specified folder and use the filename (including the `.xlsx` extension) as the genome label in the output. Annotation files should therefore be named informatively and consistently before running.

Outputs are typically renamed to reflect the function under investigation (for example, `biotin_biosynthesis_genes.csv`, `sodium-glutamate_symporter_genes.csv`) before the next search. The set of targeted output files generated this way constitutes the gene-level evidence underlying the GO enrichment results reported in the manuscript.
