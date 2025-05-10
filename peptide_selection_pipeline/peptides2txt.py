import pandas as pd
import sys

def extract_unique_peptides(input_csv, output_txt):
    try:
        df = pd.read_csv(input_csv)
        if "peptide" not in df.columns:
            print("Error: 'peptide' column not found in the CSV file.")
            sys.exit(1)

        # Extract unique, non-null, and non-empty peptides
        unique_peptides = set(df["peptide"].dropna().astype(str).str.strip())
        
        # Write unique peptides to a text file
        with open(output_txt, 'w') as f:
            for peptide in unique_peptides:
                f.write(f"{peptide}\n")
        print(f"Unique peptides written to {output_txt}")
    except FileNotFoundError:
        print(f"Error: File {input_csv} not found.")
        sys.exit(1)

def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_peptides.py <input_csv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_txt = "peptides.txt"
    extract_unique_peptides(input_csv, output_txt)

if __name__ == "__main__":
    main()
