# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 11:40:21 2020

@author: hp
"""


"""Coding peptides or 31-mers and generating data.
"""


import os
import random as rn
import numpy as np
import pandas as pd
from keras.preprocessing.sequence import pad_sequences as padding
from sklearn.model_selection import train_test_split


os.environ['PYTHONHASHSEED'] = str(0)
rn.seed(0)  # set python random seed (0)
np.random.seed(1)  # set numpy random seed (1)


# peptide / 31-mer coding
def coding(seq_type, seqs):
    # amino acid dictionaries
    if seq_type == 'peptides':
        dic = {'A':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'K':9, 
               'L':10, 'M':11, 'N':12, 'P':13, 'Q':14, 'R':15, 'S':16, 'T':17, 
               'V':18, 'W':19, 'Y':20}
    elif seq_type == '31mers':
        dic = {'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, 'H':6, 'I':7, 'K':8, 
               'L':9, 'M':10, 'N':11, 'P':12, 'Q':13, 'R':14, 'S':15, 'T':16, 
               'V':17, 'W':18, 'Y':19, 'Z':20}
    
    # coding
    coded_seqs = []
    for seq in seqs:
        coded_seq = [dic.get(aa) for aa in seq]
        coded_seqs += [coded_seq]
    
    return coded_seqs


# generate dataset
def data_encoding(pep_path='', max_len=0):
    # load peptides
    if os.path.exists(pep_path):
        print("Loading peptides.")
        df_peps = pd.read_csv(pep_path, sep='\t', names=["peptide", "label"])
        df_train, df_val = train_test_split(df_peps, test_size=0.1,
                                            random_state=2, shuffle=True,
                                            stratify=df_peps.label)
        train_peps = df_train.peptide.values
        y_train = df_train.label.values
        val_peps = df_val.peptide.values
        y_val = df_val.label.values
        
        # get the maximum length of peptides in training data set
        lengths = [len(pep) for pep in df_peps.peptide.values]
        max_len += max(lengths)
        print("The maximum length of peptides in training data is %s."
              % max_len)
        
        # data coding
        train_code = coding('peptides', train_peps)
        x_train = padding(train_code, maxlen=max_len,
                          padding='post', truncating='post', value=0)
        val_code = coding('peptides', val_peps)
        x_val = padding(val_code, maxlen=max_len, 
                        padding='post', truncating='post', value=0)
    else:
        print("The file does not exist, please check again.")
        x_train = []
        y_train = []
        x_val = []
        y_val = []
        max_len = 0
    
    return x_train, y_train, x_val, y_val, max_len