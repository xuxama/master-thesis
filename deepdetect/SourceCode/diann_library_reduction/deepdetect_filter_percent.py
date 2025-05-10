# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 15:38:42 2022

@author: hp
"""


import os
import pandas as pd
pd.set_option('display.max_columns', None)


def deepdetect_filter(fasta_file, rslt_file, regular, protease,
                      missed_cleavages, min_len, max_len,
                      speclib_file, percent):
    # fasta prediction
    print("Predicting fasta file with DeepDetect.")
    os.chdir(r"..\..\DeepDetect")
    fasta_path = f"..\SourceCode\diann_library_reduction\{fasta_file}"
    rslt_path = f"..\SourceCode\diann_library_reduction\{rslt_file}"
    cmd = (f'DeepDetect.exe --input={fasta_path} --output={rslt_path} '
           f'--regular="{regular}" --protease={protease} '
           f'--missed_cleavages={missed_cleavages} '
           f'--min_len={min_len} --max_len={max_len}')
    os.system(cmd)
    os.chdir(r"..\SourceCode\diann_library_reduction")
    
    # load predicted detectability
    print("Loading predicted detectabilities.")
    df_pred = pd.read_csv(rslt_file, sep='\t',
                          usecols=["Peptide sequence",
                                   "Peptide detectability"])
    df_pred.columns = ["PeptideSequence", "Det"]
    df_det = df_pred.groupby("PeptideSequence", as_index=False)["Det"].mean()
    
    # load spectral library
    print("Loading spectral library.")
    df_speclib = pd.read_csv(speclib_file, sep='\t')
    
    # merge data
    print("Assigining predicted results to spectral library.")
    df_merge = pd.merge(df_speclib, df_det, how='inner',
                        suffixes=('_lib', '_det'))
    if len(df_speclib) != len(df_merge):  # digestion check
        print("Warning: there are some digested peptides without predicted"
              "detectabilities!")
    df_merge["Rank"] = df_merge.Det.rank(
                            method='min', ascending=False, pct=True)
    
    # filter library
    print("Filtering spectral library.")
    df_filter = df_merge[df_merge.Rank <= percent].copy()
    df_filter.drop(["Det", "Rank"], axis=1, inplace=True)
    print("There are %s / %s = %s data in the filtered spectral library."
          % (len(df_filter), len(df_speclib),
             len(df_filter) / len(df_speclib)))
    
    # save results
    print("Saving filtered library.")
    df_filter.to_csv(r"report_lib_det%i.tsv" % (percent * 100))


if __name__ == '__main__':
    # parameters
    ## searched protein sequence file in FASTA format (.fasta)
    fasta_file = "UPS_YEAST_uniprot_sp_Release201311.fasta"
    ## output file in TXT format (.txt) of DeepDetect prediction
    rslt_file = r"deepdetect_pred.txt"
    ## the regular expression used to extract the protein id
    regular = '>(.*?)\s'
    ## the digestion protease
    protease = 'Trypsin'
    ## the maximum number of missed cleavages allowed in in silico digestion
    missed_cleavages = 2
    ## the minimum length of the theoretical peptide fragment
    min_len = 7
    ## the maximum length of the theoretical peptide fragment
    max_len = 47
    ## spectral library file in TSV format (.tsv)
    speclib_file = "report-lib.tsv"
    ## the retained percentage of peptides filtered by DeepDetect
    ## e.g., set percent=0.4 to get the top 40% spectral library
    percent = 0.4
    
    # library reduction
    deepdetect_filter(fasta_file, rslt_file, regular, protease,
                      missed_cleavages, min_len, max_len,
                      speclib_file, percent)
