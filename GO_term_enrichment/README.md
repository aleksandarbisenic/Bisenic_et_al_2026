\# Protocol for the GO term enrichment analysis



\## Generating GO annotations of bacterial genes



Bacterial genomes were annotated using PANNZER2 (http://ekhidna2.biocenter.helsinki.fi/sanspanz/) by uploading the Prokka-generated FAA files. The GO file was downloaded for further processing. The annotation file



Using the **pannzer2\_GO\_file\_parser\_script.txt** Excel Macro, the GO file was processed, so as to remove all GO terms that do not belong to the "Biological Process" ontology, producing the **parsed GO annotation file**. This was done for all of the genomes analyzed.



\## Identifying enriched GO terms in genomes of interest



For each parsed GO annotation file, the first and fourth columns were added to a new Excel file (**Fisher\_GO\_input.xlsx**), and the first of the two was named according to the name of the bacterium in question by adding a new row at the top. The second column for the first bacterium was named "GO term". Every genome occupied 2 columns.



Enrichment was calculated using the **go\_enrichment\_fisher.py** Python script. The command was as follows:



python go\_enrichment\_fisher.py Fisher\_GO\_input.xlsx NGB244 NGB245



The names of the bacteria at the end of the line define the genomes of interest, the ones whose gene enrichment the script is going to calculate compared to the rest of the genomes that exist in the **Fisher\_GO\_input.xlsx** file.



GO terms in the output file (**enriched\_GO.csv**) are ranked by decreasing p-value and FDR. Names of statistically significant GO terms were copied to a TXT list (**go\_descriptions.txt**) in order to be used for subsequent comparisons.



\## Quantification of the presence of selected genes in genomes of interest



For each genome of interest, we produced an input TXT file (e.g. NGB245\_GO\_input.txt)by copying the first four columns of the parsed GO annotation file. Python script **differential\_go\_terms.py** simply counts the number of genes assigned to GO terms listed in **go\_descriptions.txt** produced in the previous step. The script, GO terms list, and the TXT files for each genome need to be present in the same folder. The command is as follows:



python differential\_go\_terms.py go\_descriptions.txt NGB329\_GO\_input.txt NGB218\_GO\_input.txt ...



Output is called **go\_counts\_output.csv** and it's a table with gene counts per GO term per genome. All of the genomes need to be listed in the command. This kind of glutamate metabolism-related gene count comparison was used in our paper for the comparison of the genomes of selected butyrogenic gut bacteria, as well as the comparison of *Pilosibacter* genomes.

