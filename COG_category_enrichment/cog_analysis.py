import csv
from collections import defaultdict

# Predefined mapping of COG letters to functional categories
COG_CATEGORY_MAP = {
    "A": "RNA processing and modification",
    "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation, ribosomal structure, and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination, and repair",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport, and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
    "T": "Signal transduction mechanisms",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "V": "Defense mechanisms",
    "W": "Extracellular structures",
    "Y": "Nuclear structure",
    "Z": "Cytoskeleton"
}

# Input and output file paths
input_file = "cog_genes_input.txt"  # Input file with three columns: COG, Gene Description, KO Terms
enrichment_output = "cog_enrichment_results.csv"  # Enrichment results
grouped_output = "grouped_genes_by_cog.csv"  # Grouped genes by COG category

# Function to calculate COG category enrichment and group genes
def calculate_cog_enrichment_and_group(input_file):
    cog_counts = defaultdict(int)
    grouped_genes = defaultdict(list)

    # Read input file and process rows
    try:
        with open(input_file, "r") as infile:
            reader = csv.reader(infile, delimiter="\t")
            for row in reader:
                if len(row) != 3:
                    print(f"Warning: Skipping malformed row: {row}")
                    continue
                
                cog, gene_description, ko_terms = row
                for letter in cog:  # Process each letter in the COG identifier
                    if letter in COG_CATEGORY_MAP:
                        category = COG_CATEGORY_MAP[letter]
                        cog_counts[category] += 1
                        grouped_genes[category].append([letter, gene_description, ko_terms])
                    else:
                        print(f"Warning: COG identifier '{letter}' not found in category map.")
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        exit()

    return cog_counts, grouped_genes

# Calculate enrichment and group genes
cog_enrichment, grouped_genes = calculate_cog_enrichment_and_group(input_file)

# Write COG enrichment results to a CSV file
with open(enrichment_output, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["COG Category", "Description", "Gene Count"])  # Header
    for letter, description in COG_CATEGORY_MAP.items():
        count = cog_enrichment.get(description, 0)
        csvwriter.writerow([letter, description, count])

print(f"COG enrichment results saved to '{enrichment_output}'.")

# Write grouped genes by COG category to a CSV file
with open(grouped_output, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["COG Identifier", "Gene Description", "KO Terms"])  # Header
    for category, genes in grouped_genes.items():
        csvwriter.writerow([])  # Blank row to separate categories
        csvwriter.writerow([f"Category: {category}", "", ""])  # Category header
        csvwriter.writerows(genes)

print(f"Grouped genes by COG category saved to '{grouped_output}'.")
