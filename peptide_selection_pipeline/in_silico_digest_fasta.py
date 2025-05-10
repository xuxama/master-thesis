import sys
from pyteomics import parser, fasta
import pandas as pd

# Function to perform in-silico digestion of a protein sequence
def digest_protein(protein_sequence, missed_cleavages):
    return list(parser.cleave(protein_sequence, parser.expasy_rules['trypsin'], missed_cleavages=missed_cleavages))

# Function to process a FASTA file and perform digestion
def process_fasta(fasta_file, missed_cleavages, output_file):
    peptides = []
    
    with fasta.read(fasta_file) as fasta_reader:
        for entry in fasta_reader:
            protein_id, sequence = entry.description, entry.sequence
            digested_peptides = digest_protein(sequence, missed_cleavages)
            peptides.extend([(protein_id, pep) for pep in digested_peptides])
    
    # Convert to DataFrame and save to CSV
    peptides_df = pd.DataFrame(peptides, columns=['Protein_ID', 'peptide'])
    peptides_df.to_csv(output_file, index=False)
    print(f"Peptides have been saved to {output_file}.")

# Main function
def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python proteome_digest.py <FASTA file> [missed_cleavages]")
        print("  <FASTA file>      : Path to the FASTA file containing the proteome (required).")
        print("  [missed_cleavages]: Number of missed cleavages allowed (optional, default is 0).")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    missed_cleavages = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    output_file = "peptides.csv"
    
    process_fasta(fasta_file, missed_cleavages, output_file)

if __name__ == "__main__":
    main()
