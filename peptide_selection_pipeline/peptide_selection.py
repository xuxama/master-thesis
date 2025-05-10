import pandas as pd

def count_valid_peaks(peaks_str):
    """Counts the number of valid peaks from a comma-separated string."""
    if pd.isna(peaks_str) or peaks_str.strip() == "":
        return 0  # Return 0 if the string is empty or NaN
    return len(peaks_str.split(','))  # Count number of peaks

def filter_peptides(input_csv, output_csv):
    # Load the peptides.csv file
    df = pd.read_csv(input_csv)

    # Drop rows with missing values in the relevant columns
    required_columns = ['length', 'M_count', 'DMP_prob', 'DD_prob', 'CS_prob', 'valid_peaks']
    df = df.dropna(subset=required_columns)

    # Count valid peaks
    df['valid_peaks_count'] = df['valid_peaks'].apply(count_valid_peaks)

    # Base filtering conditions
    conditions = (
        (df['length'] >= 8) & (df['length'] <= 25) &
        (df['M_count'] == 0) &
        (df['DMP_prob'] >= 0.60) &
        (df['DD_prob'] >= 0.40) &
        (df['CS_prob'] >= 0.75) &
        (df['valid_peaks_count'] >= 5)
    )

    # Add unique column condition if it exists
    if 'unique' in df.columns:
        conditions = conditions & (df['unique'] != False)

    filtered_df = df[conditions]

    # Drop the temporary count column before saving
    filtered_df = filtered_df.drop(columns=['valid_peaks_count'])

    # Save the filtered dataframe to final_selection.csv
    filtered_df.to_csv(output_csv, index=False)
    print(f"Filtered peptides saved to {output_csv}")

def main():
    input_csv = "peptides.csv"
    output_csv = "final_selection.csv"
    filter_peptides(input_csv, output_csv)

if __name__ == "__main__":
    main()
