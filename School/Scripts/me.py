import sys
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
# Chemins des fichiers
gff_file = 'CDS_gene_correspondance.gff'
abundance_file = '4932-WHOLE_ORGANISM-integrated.txt'
fasta_file = 'CDS_multifasta'

def parse_gff(gff_file):
    """
    Parse the GFF file to create a dictionary mapping gene names to fasta IDs.
    """
    d_abondseq = {}

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
                    d_abondseq[gene_name] = {"fasta_id": fasta_id}
    except FileNotFoundError:
        print(f"Error: File '{gff_file}' not found.")
        sys.exit(1)
    
    return d_abondseq

def get_abundance(abundance_file, d_abondseq):
    """
    Retrieve abundance information for each CDS and store it in the dictionary.
    """
    try:
        with open(abundance_file, 'r') as file:
            for line in file:
                if line.startswith("#"):
                    continue
                data = line.strip().split("\t")
                gene_name = data[1].split(".")[1]  # Extract gene name
                abundance = data[2]
                if gene_name in d_abondseq:
                    d_abondseq[gene_name]["abundance"] = abundance
    except FileNotFoundError:
        print(f"Error: File '{abundance_file}' not found.")
        sys.exit(1)

def get_sequence(fasta_file, d_abondseq):
    try:
        with open(fasta_file, 'r') as file:
            current_seq = ""
            current_gene = ""
            for line in file:
                if line.startswith(">"):
                    if current_gene and current_seq in d_abondseq:
                        d_abondseq[current_gene]["seq"] = current_seq
                    current_gene = line.strip().split()[0][1:]
                    current_seq = ""
                else:
                    current_seq += line.strip()
            # Process the last sequence
            if current_gene and current_gene in d_abondseq:
                d_abondseq[current_gene]["seq"] = current_seq
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)

def calculate_aa_frequency(d_abondseq):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    for gene_name, info in d_abondseq.items():
        if 'seq' in info:  
            sequence = info['seq']
            length = len(sequence)
            aa_frequency = {}

            for aa in amino_acids:
                aa_frequency[aa] = 0

            for aa in sequence:
                if aa in aa_frequency:
                    aa_frequency[aa] += 1

            for aa in amino_acids:
                aa_frequency[aa] /= length

            d_abondseq[gene_name]["aa_frequency"] = aa_frequency

def plot_aa_abundance_vs_frequency(d_abondseq, amino_acid):
    x = []  # Abundance
    y = []  # Frequency

    for gene_name, info in d_abondseq.items():
        if 'aa_frequency' in info and amino_acid in info['aa_frequency'] and 'abundance' in info:
            x.append(float(info['abundance']))
            y.append(info['aa_frequency'][amino_acid])

    plt.scatter(x, y, alpha=0.5)
    plt.xlabel('Abundance')
    plt.ylabel(f'Frequency of {amino_acid}')
    plt.title(f'Abundance vs Frequency of {amino_acid}')
    plt.show()

def calculate_spearman_correlation(d_abondseq):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    for aa in amino_acids:
        x = []  # Abundance
        y = []  # Frequency
        for gene_name, info in d_abondseq.items():
            if 'aa_frequency' in info and aa in info['aa_frequency'] and 'abundance' in info:
                x.append(float(info['abundance']))
                y.append(info['aa_frequency'][aa])
        if x and y:
            correlation, p_value = spearmanr(x, y)
            print(f"Spearman correlation for {aa}: {correlation:.3f}, p-value: {p_value:.3f}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <gff_file> <abundance_file> <fasta_file>")
        sys.exit(1)

    gff_file = sys.argv[1]
    abundance_file = sys.argv[2]
    fasta_file = sys.argv[3]

    d_abondseq = parse_gff(gff_file)
    get_abundance(abundance_file, d_abondseq)
    get_sequence(fasta_file, d_abondseq)
    calculate_aa_frequency(d_abondseq)

    amino_acids_to_plot = ['E', 'D', 'K', 'R']
    for amino_acid in amino_acids_to_plot:
        plot_aa_abundance_vs_frequency(d_abondseq, amino_acid)

    calculate_spearman_correlation(d_abondseq)
