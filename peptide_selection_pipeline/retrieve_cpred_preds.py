import pandas as pd

def merge_cpred_data(peptides_csv, peptides_cpred_csv, cpred_input_csv):
    # Load dataframes
    peptides_df = pd.read_csv(peptides_csv)
    peptides_cpred_df = pd.read_csv(peptides_cpred_csv)
    cpred_input_df = pd.read_csv(cpred_input_csv)

    # Merge cpred_input with peptides_CPred in order
    merged_df = pd.concat([cpred_input_df, peptides_cpred_df], axis=1)

    # Remove 'Modifications' column and rename 'Peptide_sequence' to 'peptide'
    merged_df = merged_df.drop(columns=['Modifications'])
    merged_df = merged_df.rename(columns={'Peptide_sequence': 'peptide'})

    # Identify the column with the highest value per row (excluding 'peptide')
    probability_columns = [col for col in merged_df.columns if col.startswith('Probability_CS')]
    merged_df['CS'] = merged_df[probability_columns].idxmax(axis=1)
    merged_df['CS_prob'] = merged_df.max(axis=1, numeric_only=True)

    # Reformat 'CS' values to remove the 'Probability_CS' prefix
    merged_df['CS'] = merged_df['CS'].str.replace('Probability_CS', '')

    # Drop the original probability columns and keep relevant columns
    merged_df = merged_df[['peptide', 'CS', 'CS_prob']]

    # Merge the predictions back to the original peptides dataframe based on 'peptide'
    final_df = pd.merge(peptides_df, merged_df, on='peptide', how='left')

    # Save the updated peptides.csv
    final_df.to_csv(peptides_csv, index=False)
    print(f"Updated {peptides_csv} with CS assignments.")

def main():
    peptides_csv = "peptides.csv"
    peptides_cpred_csv = "peptides_CPred.csv"
    cpred_input_csv = "cpred_input.csv"

    merge_cpred_data(peptides_csv, peptides_cpred_csv, cpred_input_csv)

if __name__ == "__main__":
    main()