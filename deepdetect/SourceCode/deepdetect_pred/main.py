# -*- coding: utf-8 -*-
# Creation date: 20/12/2022
# Revision history:
#     20/12/2022 - Last revision
# Written by Jinghan Yang
# Copyright (c) 2022 Yan Fu's Group (http://fugroup.amss.ac.cn/)
# Academy of Mathematics and Systems Science, Chinese Academy of Sciences
# All rights reserved


"""
This program is for peptide detectability prediction for eight commonly used
proteases, i.e., Trypsin, ArgC, Chymotrypsin, GluC, LysC, AspN, LysN, and
LysargiNase. Please choose the consistent model for prediction.
"""


import time
start = time.time()
from read_fasta import read_fasta
from in_silico_digestion import digestion
from predictor import predictor
import sys
import getopt


# =============================================================================
# Get the restriction enzyme cutting site(s), cutting terminal, and the
# padding length.
# =============================================================================
def info(protease):
    if protease == 'Trypsin':
        sites = 'KR'
        terminal = 'C'
        padding_len = 51
    elif protease == 'ArgC':
        sites = 'R'
        terminal = 'C'
        padding_len = 51
    elif protease == 'Chymotrypsin':
        sites = 'WFYLM'
        terminal = 'C'
        padding_len = 51
    elif protease == 'GluC':
        sites = 'E'  # 'ED'
        terminal = 'C'
        padding_len = 48
    elif protease == 'LysC':
        sites = 'K'
        terminal = 'C'
        padding_len = 54
    elif protease == 'AspN':
        sites = 'D'
        terminal = 'N'
        padding_len = 62
    elif protease == 'LysN':
        sites = 'K'
        terminal = 'N'
        padding_len = 55
    elif protease == 'LysargiNase':
        sites = 'KR'
        terminal = 'N'
        padding_len = 43
    else:
        print("Error: This tool does not support %s protease yet,"
              " only Trypsin, ArgC, Chymotrypsin, GluC, LysC, AspN,"
              " LysN and LysargiNase are optional for now." % protease)
        sys.exit(1)
    return sites, terminal, padding_len


# =============================================================================
# Read fasta file, in silico digestion, and prediction.
# =============================================================================
def DeepDetect(data_path, res_path, regular, protease,
               missed_cleavages, min_len, max_len):
    # hyper-parameters
    sites, terminal, padding_len = info(protease)
    
    # read fasta file
    s1 = time.time()
    fasta = read_fasta(data_path, regular)
    e1 = time.time()
    print("Time cost of loading file is %s seconds." % (e1 - s1))
    
    # in silico digestion
    s2 =time.time()
    data =[]
    for protein in fasta:
        pro_id, pro_seq = protein
        digested_seqs = digestion(pro_seq, sites, terminal,
                                  missed_cleavages, min_len, max_len)
        for seqs in digested_seqs:
            data += [[pro_id, seqs]]
    e2 = time.time()
    print("Time cost of in silico digestion is %s seconds." % (e2 - s2))
    
    # prediction
    predictor(protease, data, padding_len, res_path)


if __name__ == '__main__':
    s0 = time.time()
    # default hyper-parameters
    regular = '>(.*?)\s'
    protease = 'Trypsin'
    missed_cleavages = 2
    min_len = 7
    max_len = 47
    
    # input hyper-parameters
    if len(sys.argv[1:]) <= 1:
        print("Error: Wrong command! Please read User Guide of DeepDetect.\n"
              "Example: DeepDetect.exe --input=input_filename "
              "--output=output_filename --protease='Trypsin' "
              "--missed_cleavages=2 --min_len=7 --max_len=47")
        sys.exit(1)
    else:
        options, remainder = getopt.getopt(sys.argv[1:], '',
                                           ['input=', 'output=', 'regular=',
                                           'protease=', 'missed_cleavages=',
                                           'min_len=', 'max_len='])
        for opt, arg in options:
            if opt == '--input':
                data_path = arg
            elif opt == '--output':
                res_path = arg
            elif opt == '--regular':
                regular = arg
            elif opt == '--protease':
                protease = str(arg)
            elif opt == '--missed_cleavages':
                missed_cleavages = int(arg)
            elif opt == '--min_len':
                min_len = int(arg)
            elif opt == '--max_len':
                max_len = int(arg)
            else:
                print ("Error: Argument: %s is not recognized.\n"
                       "Exiting..." % opt)
                sys.exit()
    
    # in silico digesetion and prediction
    DeepDetect(data_path, res_path, regular, protease,
               missed_cleavages, min_len, max_len)
    e0 = time.time()
    print("Time cost of the program is %s seconds." % (e0 - s0))
    
    end = time.time()
    print("Total time cost is %s seconds." % (end - start))
    print("-----The program has finished running-----")
