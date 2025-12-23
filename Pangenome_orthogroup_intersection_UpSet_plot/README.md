# Identification of shared genes and unique genes



The protocol relies on the Panaroo-generated **gene_presence_absence.csv** file. Python script **extract_panaroo.py** takes this file as input, alongside the name of one column from this file corresponding to the genome of interest, whose unique genes we wish to identify:



python extract_panaroo.py gene_presence_absence.csv NGB245



This command produces three output files: 1) **upset_matrix.tsv** (binary presence-absence table to be used in the next step), 2) **NGB245_unique_genes.tsv** (table of genes unique to NGB245, which was defined in the command as the final argument), 3) **core_genes.tsv** (list of genes that are present in all of the analyzed genomes).



Next, gene family intersections were calculated and visualized using the **plot_upset.py** Python script:



python plot_upset.py upset_matrix.tsv NGB245



The script first checks if NGB245 is present in the table (sanity check step). Then it produces an UpSet plot showing ten largest gene family intersections called **upset_pilosibacter_top10.png**.

