import pandas as pd
import sys

def format_cpred_input(peptides_csv, output_csv):
    try:
        # Load peptides.csv
        peptides_df = pd.read_csv(peptides_csv)
        
        # Keep only the 'peptide' column and rename it
        if "peptide" not in peptides_df.columns:
            print("Error: 'peptide' column not found in peptides.csv.")
            sys.exit(1)
        
        cpred_df = peptides_df[["peptide"]].rename(columns={"peptide": "Peptide_sequence"})
        
        # Print basic info before filtering
        print(f"Total peptides: {len(cpred_df)}")
        
        # Filter out peptides with length < 2
        short_peptides = cpred_df[cpred_df['Peptide_sequence'].str.len() < 2]
        if len(short_peptides) > 0:
            print(f"Removing {len(short_peptides)} peptides that are too short (< 2 characters)")
            cpred_df = cpred_df[cpred_df['Peptide_sequence'].str.len() >= 2]
            print(f"Remaining peptides after length filtering: {len(cpred_df)}")

        # Filter out peptides with non-standard amino acid "U"
        peptides_with_u = cpred_df[cpred_df['Peptide_sequence'].str.contains('U')]
        if len(peptides_with_u) > 0:
            print(f"Removing {len(peptides_with_u)} peptides containing 'U' (non-standard amino acid)")
            cpred_df = cpred_df[~cpred_df['Peptide_sequence'].str.contains('U')]

        # Add an empty 'Modifications' column
        cpred_df["Modifications"] = ""
        
        # Save to cpred_input.csv
        cpred_df.to_csv(output_csv, index=False)
        print(f"Formatted cpred_input.csv saved to {output_csv}")
    
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    peptides_csv = "peptides.csv"
    output_csv = "cpred_input.csv"
    
    format_cpred_input(peptides_csv, output_csv)

if __name__ == "__main__":
    main()
