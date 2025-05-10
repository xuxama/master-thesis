#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Filter unipept results by taxon and report problematic uniprot matches."
    )
    parser.add_argument("--taxon_id", required=True, type=int,
                        help="Taxon id to filter for (e.g. 208964)")
    parser.add_argument("--protein_id", required=True, type=str,
                        help="Protein id to exclude from problematic matches (e.g. P48632)")
    parser.add_argument("--keywords", required=True, type=str,
                        help="Comma-separated keywords (e.g. 'Ferripyoverdine receptor,Ferripyoverdine,FpvA')")
    parser.add_argument("--unipept_file", required=True, type=str,
                        help="Path to the input unipept pept2prot CSV file")
    parser.add_argument("--peptides_file", required=True, type=str,
                        help="Path to the input peptides CSV file to update")
    args = parser.parse_args()

    # Parse the comma-separated keywords and strip whitespace
    keywords = [kw.strip() for kw in args.keywords.split(",")]

    # Read the unipept results CSV
    df_uni = pd.read_csv(args.unipept_file)
    
    # Filter rows by the provided taxon_id (assumed to be in the 'taxon_id' column)
    df_uni = df_uni[df_uni["taxon_id"] == args.taxon_id]

    # Group by peptide and get unique (uniprot_id, protein_name) tuples for each peptide
    grouped = df_uni.groupby("peptide").apply(
        lambda x: x[["uniprot_id", "protein_name"]].drop_duplicates().to_dict("records")
    ).to_dict()

    # For each peptide, determine the "problematic" matches.
    # A match is considered problematic if its uniprot_id is NOT the provided protein_id
    # AND its protein_name does not contain any of the provided keywords (case-insensitive).
    problematic_matches = {}
    problematic_count = {}
    for pep, matches in grouped.items():
        problems = []
        for match in matches:
            if match["uniprot_id"] != args.protein_id:
                if not any(kw.lower() in match["protein_name"].lower() for kw in keywords):
                    problems.append(f'{match["uniprot_id"]}|{match["protein_name"]}')
        problematic_matches[pep] = problems
        problematic_count[pep] = len(problems)

    # Read peptides.csv which is assumed to have a column "peptide"
    df_pep = pd.read_csv(args.peptides_file)

    # For peptides not found in the unipept results, set as missing (pd.NA)
    def get_unique(pep):
        if pep in problematic_count:
            return problematic_count[pep] == 0
        else:
            return pd.NA

    def get_problematic(pep):
        if pep in problematic_matches:
            return ";".join(problematic_matches[pep])
        else:
            return pd.NA

    df_pep["unique"] = df_pep["peptide"].apply(get_unique)
    df_pep["problematic_uniprot_matches"] = df_pep["peptide"].apply(get_problematic)

    # Overwrite the peptides CSV file with the updated dataframe.
    df_pep.to_csv(args.peptides_file, index=False)

if __name__ == "__main__":
    main()
