import pandas as pd
import sys

def merge_peptides_with_dmp(peptides_csv, dmp_txt, output_csv):
    try:
        # Load peptides.csv
        peptides_df = pd.read_csv(peptides_csv)
        if "peptide" not in peptides_df.columns:
            print("Error: 'peptide' column not found in peptides.csv.")
            sys.exit(1)

        # Load peptides_DMP.txt (tab-separated file)
        dmp_df = pd.read_csv(dmp_txt, sep="\t")

        # Rename columns and drop the last column
        dmp_df = dmp_df.iloc[:, :2]  # Keep only the first two columns
        dmp_df.columns = ["peptide", "DMP_prob"]

        # Merge dataframes on 'peptide' allowing duplicates
        merged_df = peptides_df.merge(dmp_df, on="peptide", how="left")

        # Save updated peptides.csv
        merged_df.to_csv(output_csv, index=False)
        print(f"Updated peptides.csv saved to {output_csv}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    peptides_csv = "peptides.csv"
    dmp_txt = "peptides_DMP.txt"
    output_csv = "peptides.csv"  # Overwriting the original peptides.csv

    merge_peptides_with_dmp(peptides_csv, dmp_txt, output_csv)

if __name__ == "__main__":
    main()