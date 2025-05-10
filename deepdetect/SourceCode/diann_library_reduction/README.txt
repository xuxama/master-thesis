# Spectral library reduction by DeepDetect

This script is written in Python3.6. It can be used to reduce spectral library for DIA-NN analysis.
Specifically, the input is a tab-separated (.tsv) library file supported by DIA-NN as well as the corresponding protein sequence file (.fasta).
Then the output is the filtered library, which can be directly searched by DIA-NN with much less running time and without loss of peptide and protein identificaiton sensitivity.

## Environment requirements
python >= 3.6

tensorflow ==1.15.0 (tensorflow 1.x.x is only supported by python <= 3.7)
keras == 2.3.1
numpy == 1.19.5
pandas == 1.1.0

## Other requirement
Since the spectral library is reduced by DeepDetect, the executable program of DeepDetect must be downloaded (from http://fugroup.amss.ac.cn/software/DeepDetect/DeepDetect.html) too.
Specifically, the DeepDetect folder of the executable program and the diann_library_reduction folder of this script should be in the same path.

## Usage
1. Open the deepdetect_filter_percent.py script in the python environment.
2. Input the parameters for this program. The required parameters are: 
	fasta_file: the searched protein sequence file in FASTA format (.fasta), e.g., "UPS_YEAST_uniprot_sp_Release201311.fasta" (we used in the paper for the DIA-NN analysis)
	rslt_file: the output file in TXT format (.txt) of DeepDetect prediction, e.g., "deepdetect_pred.txt"
	regular: the regular expression used to extract the protein id, default: '>(.*?)\s'
	protease: the digestion protease, default: 'Trypsin'
	missed_cleavages: the maximum number of missed cleavages allowed in in silico digestion, default: 2
	min_len: the minimum length of the theoretical peptide fragment, default: 7
	max_len: the maximum length of the theoretical peptide fragment, default: 47
	speclib_file: the spectral library file in TSV format (.tsv), e.g., "report-lib.tsv"
	percent: the retained percentage of peptides filtered by DeepDetect, default: 0.4 (it means to get the top 40% spectral library)
3. Running the program. All the output files will be in the same folder of this script.

## Contact
If you have any questions, feedback, comments, or suggestions, please email yangjinghan@amss.ac.cn.
