import sys
from pyteomics import parser
import requests
import pandas as pd

# Function to get sequence from UniProt API
def get_protein_sequence(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        # Parse the FASTA response and extract sequence lines
        fasta_lines = response.text.splitlines()
        sequence = ''.join(line for line in fasta_lines if not line.startswith('>'))
        return sequence
    else:
        print(f"Error fetching sequence for {uniprot_id}")
        return None
# Function to digest protein sequence
def digest_protein(protein_sequence, missed_cleavages):
    peptides = list(parser.cleave(protein_sequence, parser.expasy_rules['trypsin'], missed_cleavages=missed_cleavages))
    return peptides

# Main function to orchestrate the script
def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python in_silico_digest.py <UniProt ID> [missed_cleavages]")
        print("  <UniProt ID>     : The UniProt ID of the protein to digest (required).")
        print("  [missed_cleavages]: The number of missed cleavages to allow (optional, default is 0).")
        sys.exit(1)
    
    uniprot_id = sys.argv[1]  # Get the UniProt ID from command line argument

    # If missed_cleavages argument is provided, use it, otherwise default to 0
    missed_cleavages = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    
    # Fetch the protein sequence from UniProt
    protein_sequence = get_protein_sequence(uniprot_id)
    
    # Perform in-silico digestion of the protein
    peptides = digest_protein(protein_sequence, missed_cleavages)
    
    # Create a DataFrame with a single column of peptides
    peptides_df = pd.DataFrame(peptides, columns=['peptide'])
    
    # Save the DataFrame to a CSV file
    peptides_df.to_csv('peptides.csv', index=False)
    print("Peptides have been saved to 'peptides.csv'.")

if __name__ == "__main__":
    main()