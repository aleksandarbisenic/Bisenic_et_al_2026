# Instructions for the KEGG module completeness analysis



## Generation of KO terms lists for genomes of interest



For functional annotation of the genome, eggNOG-mapper is used by uploading the FAA file produced by Prokka. Check the box marked "Proteins" and perform the annotation using the default settings.



Download the annotation Excel file. This file contains the KO terms assigned to annotated genes in the column L. Copy the entire column L and paste it into a new Excel file, then delete the first three rows (excess rows). It is in this second Excel file that the final list of KO terms belonging to the genome of a bacterium of interest will be generated.



First, using the Ctrl + F function, replace the "ko:" in all cells with "" (nothing). This deletes the "ko:". Next, copy the contents of the **CleanAndSplitColumn.txt** file. Open Visual Basic, click on the Insert option, and, in the drop menu, click on the Module option. Paste the copied macro script in here and simply close the module window and Visual basic. Execute the macro script by clicking on the Macros option and the "Run" command next to the macro name. This is the Excel macro that will clean up the KO terms list by splitting the KO terms wherever there are more than one KO term per cell. The excess KO terms are added into new cells. The final list of the selected genome's KO terms is now finalized. Copy the KO terms column (ignoring the column with newly assigned gene IDs), paste it into a TXT file, and name it appropriately. This is the input file for the KEGG module completeness Python script. 



Repeat the entire process for the rest of the genomes.



## Analysis of the completeness of KEGG modules for genomes of interest



Make sure all of the scripts and input TXT files are in the same folder. The scripts and commands were tested in an Ubuntu terminal inside Windows (WSL2). The analysis requires a stable internet connection. The first step requires the **module_completeness.py** Python script. 



Execute the following command for two example genomes - NGB244 and NGB245: 



python module_completeness.py NGB244.txt NGB245.txt



The analysis takes a while. It produces three output CSV files. First (**kegg_modules_ko_terms.csv**) simply lists the definitions of all the fetched KEGG modules and served as more of a check during the development process. Second output (**kegg_module_completeness_binary.csv**) shows the completeness of modules for the genomes in question as a binary, where 1 denotes a complete module, and 0 an incomplete one. Keep in mind that the script has a leniency rule, where a module with 3 or more steps is allowed to miss one step to still be deemed complete (purpose is to avoid false negatives due to sporadic annotation failures). Final output (**kegg_module_completeness_percentage.csv**) simply shows the percentage of steps for the module in question which the genome can fulfill through its KO terms.



Next, if the goal is to identify modules that are differentially complete in a set of genomes, we produce a completeness matrix with modules present in all of the genomes and modules present in none of the genomes removed. This requires the **keep_differentially_abundant_modules.py** Python script. Simply execute it in the same folder as the previous three CSV outputs:



python keep_differentially_abundant_modules.py



This produces two completeness outputs of the same structure as the previous ones, but with fewer modules, only the ones that are differentially complete. 



To produce the final heatmap, simply execute the following command in the same folder as the outputs of the previous step using the **heatmap_binary.py** Python script:



python heatmap_binary.py



The heatmap generating script allows for customization in terms of color mapping and graph dimensions.



