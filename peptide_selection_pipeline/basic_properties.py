import pandas as pd
import sys

def process_peptides():
    # Load the peptides.csv file
    try:
        peptides_df = pd.read_csv("peptides.csv")
    except FileNotFoundError:
        print("Error: peptides.csv not found. Make sure you have run in_silico_digest.py first.")
        sys.exit(1)

    # Adding new columns
    peptides_df["length"] = peptides_df["peptide"].apply(len)
    peptides_df["M_count"] = peptides_df["peptide"].str.count('M')

    # Save the updated DataFrame back to peptides.csv
    peptides_df.to_csv("peptides.csv", index=False)

    print("Updated peptides.csv with new columns: 'length' and 'M_count'.")

if __name__ == "__main__":
    process_peptides()