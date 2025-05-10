import pandas as pd
import numpy as np

def load_spectronaut_tsv(tsv_file):
    """
    Load the Spectronaut TSV file and remove rows with zero relative intensity.
    """
    df = pd.read_csv(tsv_file, sep="\t")
    df = df[df['RelativeFragmentIntensity'] > 0]
    return df

def find_interfered_peaks(peptide_df, mz_separation_threshold):
    """
    Mark peaks in the peptide group that are interfered by other peaks within the m/z threshold.
    Returns a set of indices of interfered peaks.
    """
    mzs = peptide_df['FragmentMz'].values
    interfered = set()

    for i in range(len(mzs)):
        for j in range(len(mzs)):
            if i == j:
                continue
            if abs(mzs[i] - mzs[j]) < mz_separation_threshold:
                interfered.add(i)
                break  # no need to check further once interfered

    return interfered

def process_spectronaut(tsv_file, intensity_threshold=0.05, mz_separation_threshold=10, top_n=6):
    """
    Process the Spectronaut file to compute valid y-ion peaks and average intensity for each peptide.
    Returns a DataFrame with columns: StrippedPeptide, ValidPeaks, AvgIntensity.
    """
    df = load_spectronaut_tsv(tsv_file)
    df['peak_id'] = df['FragmentType'] + df['FragmentNumber'].astype(int).astype(str)

    results = []

    for peptide, group in df.groupby('StrippedPeptide'):
        group = group.reset_index(drop=True)

        # Detect interfered peaks using all fragment ions
        interfered_indices = find_interfered_peaks(group, mz_separation_threshold)

        # Mark interference
        group['is_interfered'] = group.index.isin(interfered_indices)

        # Filter to valid peaks: y-ions, FragmentNumber > 2, above intensity threshold, not interfered
        valid_peaks = group[
            (group['FragmentType'] == 'y') &
            (group['FragmentNumber'] > 2) &
            (group['RelativeFragmentIntensity'] > intensity_threshold) &
            (~group['is_interfered'])
        ]

        valid_peak_list = valid_peaks['peak_id'].tolist()
        num_valid_peaks = len(valid_peaks)

        if num_valid_peaks >= 5:
            top_peaks = valid_peaks.nlargest(top_n, 'RelativeFragmentIntensity')
            avg_intensity = top_peaks['RelativeFragmentIntensity'].mean()
        else:
            avg_intensity = np.nan

        # Join peaks with commas, wrap the entire result in quotes for proper CSV formatting
        valid_peaks_str = '"{}"'.format(', '.join(valid_peak_list))

        results.append({
            'StrippedPeptide': peptide,
            'ValidPeaks': valid_peaks_str,
            'AvgIntensity': avg_intensity
        })

    return pd.DataFrame(results)

def main():
    peptides_csv = "peptides.csv"
    spectronaut_tsv = "peptides_ms2pip.spectronaut.tsv"

    peptides_df = pd.read_csv(peptides_csv)
    spec_results = process_spectronaut(spectronaut_tsv)

    merged_df = pd.merge(peptides_df, spec_results,
                         left_on='peptide',
                         right_on='StrippedPeptide',
                         how='left')

    merged_df = merged_df.drop(columns=['StrippedPeptide'])
    merged_df = merged_df.rename(columns={
        'ValidPeaks': 'valid_peaks',
        'AvgIntensity': 'average_peak_intensity'
    })

    merged_df.to_csv(peptides_csv, index=False)
    print(f"Updated {peptides_csv} with 'valid_peaks' and 'average_peak_intensity' columns.")

if __name__ == "__main__":
    main()
