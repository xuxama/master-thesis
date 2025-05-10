import pandas as pd

def generate_ms2pip_input(peptides_csv, output_tsv):
    # Load the peptides.csv file
    peptides_df = pd.read_csv(peptides_csv)
    
    # Extract 'peptide' and 'CS' columns
    if 'peptide' not in peptides_df.columns or 'CS' not in peptides_df.columns:
        raise ValueError("peptides.csv must contain 'peptide' and 'CS' columns.")
    
    # Drop rows where 'CS' is empty
    peptides_df = peptides_df.dropna(subset=['CS'])
    
    # Ensure 'CS' is formatted correctly by removing '+' and converting to integer
    peptides_df['CS'] = peptides_df['CS'].astype(str).str.replace('+', '').astype(float).astype(int).astype(str)
    peptides_df['peptidoform'] = peptides_df['peptide'] + '/' + peptides_df['CS']
    
    # Add spectrum_id as a sequential number starting from 1
    peptides_df['spectrum_id'] = range(1, len(peptides_df) + 1)
    
    # Select relevant columns
    ms2pip_df = peptides_df[['peptidoform', 'spectrum_id']]
    
    # Save to a TSV file
    ms2pip_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Generated {output_tsv} for MS2PIP input.")

def main():
    peptides_csv = "peptides.csv"
    output_tsv = "ms2pip_input.tsv"
    generate_ms2pip_input(peptides_csv, output_tsv)

if __name__ == "__main__":
    main()
