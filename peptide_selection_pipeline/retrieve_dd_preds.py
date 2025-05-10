import pandas as pd
import sys

def merge_dd_prob(peptides_csv, dd_txt):
    try:
        # Load peptides.csv
        peptides_df = pd.read_csv(peptides_csv)
        if "peptide" not in peptides_df.columns:
            print("Error: 'peptide' column not found in peptides.csv.")
            sys.exit(1)

        # Load peptides_DD.txt (tab-separated file)
        dd_df = pd.read_csv(dd_txt, sep="\t")
        required_dd_cols = {"Protein id", "Peptide sequence", "Peptide detectability"}
        if not required_dd_cols.issubset(dd_df.columns):
            print("Error: Required columns not found in peptides_DD.txt.")
            sys.exit(1)

        # Rename columns in dd_df to standard names
        dd_df = dd_df.rename(columns={
            "Protein id": "Protein_ID",
            "Peptide sequence": "peptide",
            "Peptide detectability": "DD_prob"
        })

        # Determine merge strategy based on peptides_df columns
        if "Protein_ID" in peptides_df.columns:
            merge_cols = ["Protein_ID", "peptide"]
        else:
            merge_cols = ["peptide"]
            # If Protein_ID is not present in peptides_df, drop it from dd_df to avoid merging it in.
            dd_df = dd_df.drop("Protein_ID", axis=1)

        # Merge dataframes on the determined columns
        merged_df = peptides_df.merge(dd_df, on=merge_cols, how="left")

        # Overwrite peptides.csv with the updated data
        merged_df.to_csv(peptides_csv, index=False)
        print(f"Updated {peptides_csv} with DD_prob values.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    peptides_csv = "peptides.csv"
    dd_txt = "peptides_DD.txt"
    merge_dd_prob(peptides_csv, dd_txt)

if __name__ == "__main__":
    main()
