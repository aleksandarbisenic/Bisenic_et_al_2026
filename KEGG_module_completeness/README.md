# KEGG module completeness analysis

This subdirectory contains the scripts and protocol for the KEGG module completeness analysis described in §2.10.1 of Bisenić et al. (2026), corresponding to Figure 2B.

## Software requirements

- Python ≥3.8
- Python packages: `pandas`, `requests`, `numpy`, `scipy`, `statsmodels`
- Microsoft Excel with VBA enabled, for KO term list preparation
- Stable internet connection, for KEGG REST API queries in step 2
- GraphPad Prism v10, for the final heatmap

Scripts and commands were tested on Ubuntu inside Windows (WSL2).

## Files in this directory

- `module_completeness.py` — fetches KEGG modules from the REST API and computes per-genome completeness
- `keep_differentially_abundant_modules.py` — filters to modules that differ across genomes
- `fisher_module_completeness.py` — Fisher's exact test for differential completeness between two groups of genomes
- `CleanAndSplitColumn.txt` — Excel VBA macro for preparing the KO term input lists

## Workflow

### Step 1. Prepare KO term lists from eggNOG-mapper annotations

For each genome of interest:

1. Submit the Prokka-generated `.faa` protein file to eggNOG-mapper, check the "Proteins" box, and run the annotation with default settings.
2. Download the annotation Excel file. KO terms are listed in column **L** (`KEGG_ko`).
3. In a new Excel workbook, paste column L into column A. Delete the first three header rows.
4. With column A selected, open Find & Replace (Ctrl + H). Enter `ko:` in "Find what", leave "Replace with" empty, and click "Replace All".
5. Open the Visual Basic editor (Alt + F11). From the menu, select Insert > Module. Paste the entire contents of `CleanAndSplitColumn.txt` into the module window, then close the editor.
6. Run the macro by pressing Alt + F8, selecting `CleanSplitAndNumberGenes`, and clicking Run. The macro splits comma- and slash-separated multi-KO entries into individual rows, removes empty rows, places one KO term per row in column B, and adds sequential gene labels in column A.
7. Copy column B only (ignoring the gene labels in column A) into a plain text file named after the genome (e.g., `NGB245.txt`). The file should contain one KO term per line, with no header.

Repeat for each genome to be compared.

### Step 2. Compute module completeness

Place all KO term `.txt` files and the Python scripts in the same working directory and run:

```
python module_completeness.py NGB245.txt P_fragilis.txt A_hadrus.txt [...]
```

Pass each genome's KO term file as a separate argument. The script fetches all KEGG modules (M00001–M01000) via the KEGG REST API, parses their step definitions, evaluates each step against each genome's KO list, and writes three output CSVs:

- `kegg_module_ko_terms.csv` — fetched module names and their step definitions. This file was produced as a sanity check during development and is not required downstream.
- `kegg_module_completeness_binary.csv` — binary matrix (rows = genomes, columns = modules), where `1` denotes a complete module and `0` an incomplete one. **Leniency rule**: modules with three or more steps are deemed complete if at most one step is missing, to mitigate sporadic annotation failures.
- `kegg_module_completeness_percentage.csv` — percentage of fulfilled steps per genome × module.

This step is the slowest part of the analysis (approximately 30–60 minutes depending on network conditions) and requires uninterrupted internet access for the duration of the run.

### Step 3. Filter to differentially-present modules

Run in the same directory as the step 2 outputs:

```
python keep_differentially_abundant_modules.py
```

The script reads the binary and percentage matrices, removes modules that are complete in all genomes or incomplete in all genomes (which carry no comparative signal), and — when more than 10 genomes are included — additionally removes modules that differ in only a single genome (noise filter). It writes:

- `differentially_present_binary.csv`
- `differentially_present_percentage.csv`

### Step 4. Test for differential completeness between two groups

Run Fisher's exact test on the filtered binary matrix, specifying the two groups of genomes to compare. For the manuscript, group A comprised the saccharolytic butyrogen reference genomes and group B comprised the *Pilosibacter* genomes:

```
python fisher_module_completeness.py \
  --in differentially_present_binary.csv \
  --out fisher_results.csv \
  --groupA A_hadrus A_hallii A_rectalis C_catus C_eutactus F_prausnitzii R_hominis R_intestinalis R_inulinivorans S_variabile \
  --groupB NGB245 P_fragilis P_sp900066055 P_sp902496665 P_sp934370915 P_sp963622995 P_sp963646925
```

Genome names passed to `--groupA` and `--groupB` must match the file basenames (without the `.txt` extension) used in step 2. Default settings apply a two-sided Fisher's exact test followed by Benjamini–Hochberg correction across all tested modules. The output `fisher_results.csv` contains, per module: completeness counts in each group, group prevalences, the prevalence difference (group A − group B), the odds ratio, the raw p-value, and the BH q-value.

Genome accessions used for group A and group B in this study are listed in Supplementary Table 2 of Bisenić et al. (2026).

### Step 5. Heatmap (Figure 2B)

Modules with q ≤ 0.05 in `fisher_results.csv` were used to select the corresponding columns of `kegg_module_completeness_binary.csv`. The resulting binary matrix was imported into GraphPad Prism v10 to generate the heatmap shown in Figure 2B.
