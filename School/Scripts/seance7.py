import sys
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

d_abondseq = {}
## Créez une fonction (Parse_gff()) qui lira le fichier « CDS_gene_correspondance.gff ». 
def parse_gff(gff_file):
    """
    Parse the GFF file to create a dictionary mapping gene names to fasta IDs.
    """
    gff_file = "CDS_gene_correspondance.gff"

    try:
        with open(gff_file, 'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue
                data = line.strip().split("\t")
                gene_name = None
                fasta_id = None
                for field in data[8].split(";"):
                    key, value = field.split("=")
                    if key == "Name":
                        gene_name = value
                    elif key == "ID":
                        fasta_id = value
                if gene_name and fasta_id:
                    d_abondseq[gene_name] = fasta_id
    except FileNotFoundError:
        print(f"Error: File '{gff_file}' not found.")
        sys.exit(1)
    
    return d_abondseq
print(d_abondseq)
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <gff_file> <abundance_file> <fasta_file>")
        sys.exit(1)

    gff_file = sys.argv[1]
    abundance_file = sys.argv[2]
    fasta_file = sys.argv[3]

    d_abondseq = parse_gff(gff_file)  # Appel de la fonction parse_gff

   