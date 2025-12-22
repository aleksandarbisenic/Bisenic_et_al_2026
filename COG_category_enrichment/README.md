# Protocol for COG functional category enrichment



## Input file preparation



Input files for this analysis were prepared from the eggNOG-mapper annotation outputs for each bacterial genome separately. From the annotation file, columns G, H, and L were copied and pasted to a new Excel file in this order. Then the first three rows were deleted (they contain unnecessary headers). Next, using Ctrl + F, "ko:" was replaced with "" (empty space) in order to clean up the columns. Finally, the three columns generated in this way were copied and pasted to the final TXT input file called **cog_genes_input.txt**. All files for this and the next step need to exist in the same folder. 



## COG category gene classification



The table summarizing gene counts for each of the COG functional categories was also produced for each genome separately. Python script **cog_analysis.py** was executed in the same folder as **cog_genes_input.txt** generated for the genome in question:



python cog_analysis.py 



This script produces two outputs: 1) **cog_enrichment_results.csv**, which contain gene counts for each COG category, and 2) **grouped_genes_by_cog.csv**, which only serve for additional overview of the genes that had been assigned to each COG category.



These first two steps were performed for all six isolates analyzed in our study.



## COG category enrichment using Fisher's over-representation test



A summary Excel table was produced from **cog_enrichment_results.csv** files of all of the genomes to be analyzed and named **cog_enrichment_input.xlsx**. First column contained the names of all of the COG categories in the same order as in the **cog_enrichment_results.csv** file. Subsequent columns contained gene counts for each of the genomes copied from their respective **cog_enrichment_results.csv** files. An additional row was added at the top for column names. First column was named "COG Category" and the rest of the columns were named according to the preferred names for the analyzed genomes. Then, the enrichment was calculated for genomes of interest compared to the rest of the genomes in **cog_enrichment_input.xlsx** using the following command:



python cog_enrichment_fisher.py cog_enrichment_input.xlsx NGB244 NGB245



The command identifies COG categories that are enriched in the genomes of NGB244 and NGB245 compared to the remaining four genomes. COG categories are ranked by decreasing p-value and FDR. Output file is **cog_enrichment_fisher.csv**.

