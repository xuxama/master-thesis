# Loading modules
import pandas as pd
import numpy as np
import regex as re
from tensorflow import keras


# Dictionaries
Isoelectric_point = {
    "A": 6.00,
    "C": 5.07,
    "D": 2.77,
    "E": 3.22,
    "F": 5.48,
    "G": 5.97,
    "H": 7.59,
    "I": 6.02,
    "K": 9.74,
    "L": 5.98,
    "M": 5.74,
    "N": 5.41,
    "P": 6.30,
    "Q": 5.65,
    "R": 10.76,
    "S": 5.68,
    "T": 5.60,
    "V": 5.96,
    "W": 5.89,
    "Y": 5.66,
    "O": 9.74
}

Hydrophobicity = {      # Black et al. 1991 https://doi.org/10.1016/0003-2697(91)90045-U
    "A": 0.616,
    "C": 0.680,
    "D": 0.028,
    "E": 0.043,
    "F": 1.000,
    "G": 0.501,
    "H": 0.165,
    "I": 0.943,
    "K": 0.283,
    "L": 0.943,
    "M": 0.738,
    "N": 0.236,
    "P": 0.711,
    "Q": 0.251,
    "R": 0.00001,
    "S": 0.359,
    "T": -0.70,
    "V": 0.825,
    "W": 0.878,
    "Y": 0.880,
    "O": 0.283
}

Sequence_One_Hot = {'A': 1,
                    'C': 2,
                    'D': 3,
                    'E': 4,
                    'F': 5,
                    'G': 6,
                    'H': 7,
                    'I': 8,
                    'K': 9,
                    'L': 10,
                    'M': 11,
                    'N': 12,
                    'P': 13,
                    'Q': 14,
                    'R': 15,
                    'S': 16,
                    'T': 17,
                    'V': 18,
                    'W': 19,
                    'Y': 20,
                    'O': 21}

Charge_State_OneHot = {
    1: [1, 0, 0, 0, 0, 0],
    2: [0, 1, 0, 0, 0, 0],
    3: [0, 0, 1, 0, 0, 0],
    4: [0, 0, 0, 1, 0, 0],
    5: [0, 0, 0, 0, 1, 0],
    6: [0, 0, 0, 0, 0, 1]
}

ElementalComp = {
    'Amino Acid': ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'O'],
    'C': [3, 6, 4, 4, 3, 5, 5, 2, 6, 6, 6, 6, 5, 9, 5, 3, 4, 11, 9, 5, 12],
    'H': [5, 12, 6, 5, 5, 7, 8, 3, 7, 11, 11, 12, 9, 9, 7, 5, 7, 10, 9, 9, 21],
    'N': [1, 4, 2, 1, 1, 1, 2, 1, 3, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 3],
    'O': [1, 1, 2, 3, 1, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 3],
    'S': [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
}
ElementalComp = pd.DataFrame(ElementalComp)
ElementalComp = ElementalComp.sort_values(by='Amino Acid')
ElementalComp.set_index('Amino Acid', inplace=True)
ElementalComp.index.name = 'Amino Acid'


# All functions related to feature engineering
def tryptic(dataframe):
    dataframe["Tryptic"] = list(map(lambda x: x.endswith('R') | x.endswith('K'), dataframe['Peptide_sequence']))
    return dataframe


def peptide_length(dataframe):
    dataframe["Peptide_Length"] = dataframe["Peptide_sequence"].str.len()
    return dataframe


def basic(dataframe):
    dataframe["Fraction_Basic"] = dataframe['Peptide_sequence'].str.count('[RHKO]') / dataframe['Peptide_Length']
    return dataframe


def acid(dataframe):
    dataframe["Fraction_Acidic"] = dataframe['Peptide_sequence'].str.count('[DE]') / dataframe['Peptide_Length']
    return dataframe


def non_polar(dataframe):
    dataframe["Fraction_Non_polar"] = dataframe['Peptide_sequence'].str.count('[AFGILMPVW]') / dataframe['Peptide_Length']
    return dataframe


def polar(dataframe):
    dataframe["Fraction_Polar"] = dataframe['Peptide_sequence'].str.count('[CNQSTY]') / dataframe['Peptide_Length']
    return dataframe


def polar_basic(dataframe):
    dataframe["Fraction_Polar_basic"] = dataframe['Peptide_sequence'].str.count('[HKRO]') / dataframe['Peptide_Length']
    return dataframe


def polar_acidic(dataframe):
    dataframe["Fraction_Polar_acidic"] = dataframe['Peptide_sequence'].str.count('[DE]') / dataframe['Peptide_Length']
    return dataframe


def aromatic(dataframe):
    dataframe["Fraction_Aromatic"] = dataframe['Peptide_sequence'].str.count('[FWY]') / dataframe['Peptide_Length']
    return dataframe


def alpha(dataframe):
    dataframe["Fraction_Alpha"] = dataframe['Peptide_sequence'].str.count('[FILVWY]') / dataframe['Peptide_Length']
    return dataframe


def beta(dataframe):
    dataframe["Fraction_Beta"] = dataframe['Peptide_sequence'].str.count('[GNPS]') / dataframe['Peptide_Length']
    return dataframe


def turn(dataframe):
    dataframe["Fraction_Turn"] = dataframe['Peptide_sequence'].str.count('[AELM]') / dataframe['Peptide_Length']
    return dataframe


def sequence_checker(sequence):
    if sequence.endswith('R') | sequence.endswith('K'):
        Arginines = len(re.findall(pattern='R', string=sequence))
        Lysines = len(re.findall(pattern='K', string=sequence))
        Inhibited_RP = len(re.findall(pattern='RP', string=sequence))
        Inhibited_KP = len(re.findall(pattern='KP', string=sequence))
        Missed_Cleavage = Arginines + Lysines - Inhibited_RP - Inhibited_KP - 1
        return Missed_Cleavage
    else:
        Missed_Cleavage = 0
        return Missed_Cleavage


def missed_cleavages(dataframe):
    dataframe["Missed_Cleavages"] = dataframe['Peptide_sequence'].apply(lambda x: sequence_checker(x))
    return dataframe


def isoelectric_list(sequence):
    storage = []
    for AA in sequence:
        isoelectric_value = Isoelectric_point[AA]
        storage.append(isoelectric_value)
    return storage


def isoelectric_point(dataframe):
    dataframe["Isoelectric_point_AA"] = dataframe['Peptide_sequence'].apply(lambda x: isoelectric_list(x))
    return dataframe


def average_isoelectric_point(dataframe):
    dataframe["Avg_Isoelectric_point_AA"] = dataframe['Isoelectric_point_AA'].apply(lambda x: np.mean(x))
    return dataframe


def hydrophobicity_list(sequence):
    storage = []
    for AA in sequence:
        Hydrophobicity_value = Hydrophobicity[AA]
        storage.append(Hydrophobicity_value)
    return storage


def hydrophobicity(dataframe):
    dataframe["Hydrophobicity_AA"] = dataframe['Peptide_sequence'].apply(lambda x: hydrophobicity_list(x))
    return dataframe


def average_hydrophobicity(dataframe):
    dataframe["Avg_Hydrophobicity_AA"] = dataframe['Hydrophobicity_AA'].apply(lambda x: np.mean(x))
    return dataframe


# Define all unique elements as a global variable
Unique_elements = ['N', 'HexN', 'C', 'B', 'Li', 'Se', 'K', 'Zn', 'NeuGc', 'HexA', 'Ag', 'Hep', 'Sulf', '18O', 'Fe', '15N',
                   '2H', 'S', 'I', 'Hex', 'Ac', 'Na', 'O', 'Hg', 'F', 'Ni', 'As', 'Br', 'Mo', 'Me', 'Kdn', 'Ca', 'dHex', 'Pent',
                   'Al', 'HexNAc', 'Cu', 'H', 'Cl', 'P', '13C', 'Mg', 'NeuAc']

def elemental_comp(dataframe):
    print(f"Starting elemental composition calculation for {len(dataframe)} peptides...")
    
    # Initialize all elements in Unique_elements with 0
    for column in Unique_elements:
        dataframe[column] = 0
    
    # Create element dictionary for fast lookup (only for elements that exist in ElementalComp)
    elements_to_calculate = ["C", "H", "O", "N", "S"]  # Elements that we can actually calculate
    
    print("Creating element lookup dictionary...")
    element_dict = {}
    for amino_acid in ElementalComp.index:
        element_dict[amino_acid] = {}
        for element in elements_to_calculate:
            element_dict[amino_acid][element] = ElementalComp.loc[amino_acid][element]
    
    print("Processing peptide elemental compositions...")
    
    # Process in batches
    batch_size = max(1, len(dataframe) // 20)
    total_batches = (len(dataframe) + batch_size - 1) // batch_size
    
    for batch_num, i in enumerate(range(0, len(dataframe), batch_size)):
        batch_end = min(i + batch_size, len(dataframe))
        print(f"Processing batch {batch_num+1}/{total_batches} (peptides {i+1}-{batch_end})...")
        
        for idx, row in dataframe.iloc[i:batch_end].iterrows():
            peptide_sequence = row['Peptide_sequence']
            
            # Calculate elements for this peptide (only those in ElementalComp)
            for element in elements_to_calculate:
                element_count = 0
                
                # Sum elements from each amino acid
                for amino_acid in peptide_sequence:
                    if amino_acid in element_dict:
                        element_count += element_dict[amino_acid][element]
                
                # Add water molecule
                if element == "H":
                    element_count += 2
                elif element == "O":
                    element_count += 1
                    
                # Update the dataframe
                dataframe.at[idx, element] = element_count
    
    print("Completed elemental composition calculation")
    return dataframe

# Modified functions elemental_comp_modifications and apply_modifications to improve efficiency for large datasets 
# Avoids redownloading the file for each row and uses more efficient lookup methods

def elemental_comp_modifications(dataframe):
    print(f"Starting elemental_comp_modifications for {len(dataframe)} rows")
    
    # Download modifications file once outside the apply function
    try:
        print("Downloading modifications file...")
        
        # Load modifications data once
        modifications_unimod = pd.read_excel(
            "https://raw.githubusercontent.com/VilenneFrederique/CPred/master/CPred/Data/Unimod_modifications.xlsx")
        
        print(f"Downloaded {len(modifications_unimod)} modification entries")
        
        # Create a dictionary for faster lookups
        modifications_dict = {}
        for _, row in modifications_unimod.iterrows():
            modifications_dict[row['name']] = row
        
        print(f"Created lookup dictionary with {len(modifications_dict)} entries")
        
        # Define a modified apply_modifications function that uses pre-loaded data
        def apply_modifications_optimized(row):
            nonlocal processed_rows
            
            # Print progress every 500 rows
            if processed_rows % 500 == 0:
                print(f"Processing row {processed_rows}/{len(dataframe)}")
            
            modifications_string = row['Modifications']
            if pd.isna(modifications_string):
                processed_rows += 1
                return row
            else:
                modifications_list = re.findall(pattern="\|([\da-zA-Z]+)", string=modifications_string)
                
                # Occasionally print some details about modifications found
                if processed_rows % 1000 == 0 and modifications_list:
                    print(f"Sample modifications in row {processed_rows}: {modifications_list[:3]}")
                
                for modification in modifications_list:
                    if modification in modifications_dict:
                        modification_row = modifications_dict[modification]
                        for element in Unique_elements:
                            row[element] += modification_row[element]
                
                processed_rows += 1
                return row
        
        # Apply modifications with the pre-loaded data
        processed_rows = 0
        print("Starting to apply modifications...")
        
        # Check if Unique_elements is defined
        if 'Unique_elements' not in globals():
            print("Warning: Unique_elements not defined!")
        
        dataframe_modified = dataframe.apply(apply_modifications_optimized, axis=1)
        
        print("Finished applying modifications")
        return dataframe_modified
    
    except Exception as e:
        print(f"Error loading modifications: {e}")
        print("Continuing without modifications...")
        return dataframe


def N_term_AA_Checker(sequence):
    N_term_AA_sequence = sequence[0]
    One_hot_N_Term_AA = Sequence_One_Hot[N_term_AA_sequence]
    return One_hot_N_Term_AA


def N_term_AA(dataframe):
    dataframe["N_Term_AA"] = dataframe['Peptide_sequence'].apply(lambda x: N_term_AA_Checker(x))
    return dataframe


def C_term_AA_Checker(sequence):
    C_term_AA_sequence = sequence[-1]
    One_hot_C_Term_AA = Sequence_One_Hot[C_term_AA_sequence]
    return One_hot_C_Term_AA


def C_term_AA(dataframe):
    dataframe["C_Term_AA"] = dataframe['Peptide_sequence'].apply(lambda x: C_term_AA_Checker(x))
    return dataframe


def N_term_AA_pI_checker(sequence):
    N_term_AA_sequence = sequence[0]
    pI_N_Term_AA = Isoelectric_point[N_term_AA_sequence]
    return pI_N_Term_AA


def N_term_AA_pI(dataframe):
    dataframe["N_Term_AA_pI"] = dataframe['Peptide_sequence'].apply(lambda x: N_term_AA_pI_checker(x))
    return dataframe


def N_term_AA_pI_checker_last_2(sequence):
    N_term_AA_sequence_AA_0 = sequence[0]
    N_term_AA_sequence_AA_1 = sequence[1]
    pI_N_Term_AA_0 = Isoelectric_point[N_term_AA_sequence_AA_0]
    pI_N_Term_AA_1 = Isoelectric_point[N_term_AA_sequence_AA_1]
    pI_N_AVG = (pI_N_Term_AA_0 + pI_N_Term_AA_1) / 2
    return pI_N_AVG


def N_term_AA_pI_last_2(dataframe):
    dataframe["N_Term_AA_pI_last_2_avg"] = dataframe['Peptide_sequence'].apply(lambda x: N_term_AA_pI_checker_last_2(x))
    return dataframe


def C_term_AA_pI_checker(sequence):
    C_term_AA_sequence = sequence[-1]
    pI_C_Term_AA = Isoelectric_point[C_term_AA_sequence]
    return pI_C_Term_AA


def C_term_AA_pI(dataframe):
    dataframe["C_Term_AA_pI"] = dataframe['Peptide_sequence'].apply(lambda x: C_term_AA_pI_checker(x))
    return dataframe


def C_term_AA_pI_checker_last_2(sequence):
    C_term_AA_sequence_AA_0 = sequence[-1]
    C_term_AA_sequence_AA_1 = sequence[-2]
    pI_C_Term_AA_0 = Isoelectric_point[C_term_AA_sequence_AA_0]
    pI_C_Term_AA_1 = Isoelectric_point[C_term_AA_sequence_AA_1]
    pI_C_AVG = (pI_C_Term_AA_0 + pI_C_Term_AA_1) / 2
    return pI_C_AVG


def C_term_AA_pI_last_2(dataframe):
    dataframe["C_Term_AA_pI_last_2_avg"] = dataframe['Peptide_sequence'].apply(lambda x: C_term_AA_pI_checker_last_2(x))
    return dataframe


def monoisotopic_mass(dataframe):
    # Define mass dictionary
    monoisotopic_masses = {
        "N": 14.003074, "HexN": 161.068808, "C": 12, "B": 11.0093055, "Li": 7.016003,
        "Se": 79.9165196, "K": 38.9637074, "Zn": 63.9291448, "NeuGc": 307.090331, "HexA": 176.032088,
        "Ag": 106.905092, "Hep": 192.063388, "Sulf": 31.9720707, "18O": 17.9991603, "Fe": 55.9349393,
        "15N": 15.00010897, "2H": 2.014101779, "S": 31.9720707, "I": 126.904473, "Hex": 162.052824,
        "Ac": 59.013851, "Na": 22.9897677, "O": 15.99491463, "Hg": 201.970617, "F": 18.99840322,
        "Ni": 57.9353462, "As": 74.9215942, "Br": 78.9183361, "Mo": 97.9054073, "Me": 15.023475105,
        "Kdn": 268.079437, "Ca": 39.9625906, "dHex": 146.057909, "Pent": 132.042259, "Al": 26.9815386,
        "HexNAc": 203.079373, "Cu": 62.9295989, "H": 1.007825035, "Cl": 34.96885272, "P": 30.973762,
        "13C": 13.00335483, "Mg": 23.9850423, "NeuAc": 291.095417
    }
    
    # Use vectorized operations instead of loops
    # Initialize with zeros
    dataframe["Monoisotopic_mass"] = 0.0
    
    # Add each component's contribution
    for element, mass in monoisotopic_masses.items():
        if element in dataframe.columns:
            dataframe["Monoisotopic_mass"] += dataframe[element] * mass
    
    return dataframe


def average_mass(dataframe):
    # Define mass dictionary
    average_masses = {
        "N": 14.0067, "HexN": 161.1558, "C": 12.0107, "B": 10.811, "Li": 6.941,
        "Se": 78.96, "K": 39.0983, "Zn": 65.409, "NeuGc": 307.2540, "HexA": 176.1241,
        "Ag": 107.8682, "Hep": 192.1666, "Sulf": 32.065, "18O": 17.9991603, "Fe": 55.845,
        "15N": 15.00010897, "2H": 2.014101779, "S": 32.065, "I": 126.90447, "Hex": 162.1406,
        "Ac": 59.045, "Na": 22.98977, "O": 15.9994, "Hg": 200.59, "F": 18.9984032,
        "Ni": 58.6934, "As": 74.9215942, "Br": 79.904, "Mo": 95.94, "Me": 15.03452,
        "Kdn": 268.218, "Ca": 40.078, "dHex": 146.1412, "Pent": 132.1146, "Al": 26.9815386,
        "HexNAc": 203.1925, "Cu": 63.546, "H": 1.00794, "Cl": 35.453, "P": 30.973761,
        "13C": 13.00335483, "Mg": 24.305, "NeuAc": 291.2546
    }
    
    # Use vectorized operations instead of loops
    # Initialize with zeros
    dataframe["Average_mass"] = 0.0
    
    # Add each component's contribution
    for element, mass in average_masses.items():
        if element in dataframe.columns:
            dataframe["Average_mass"] += dataframe[element] * mass
    
    return dataframe


def encoder_sequence(sequence):
    encoded_sequence = [Sequence_One_Hot[aa] for aa in sequence]
    return encoded_sequence


def sequence_onehot(dataframe):
    # First convert sequences to one-hot encoding
    dataframe["Sequences_OneHot"] = dataframe.apply(lambda x: encoder_sequence(x["Peptide_sequence"]), axis=1)
    
    # Pad sequences
    Encoded_sequences = dataframe["Sequences_OneHot"].to_numpy()
    Encoded_sequences = keras.preprocessing.sequence.pad_sequences(
        Encoded_sequences,
        maxlen=50,
        padding='post',
        truncating='post'
    )
    
    # Store the padded sequences
    dataframe["Sequences_OneHot"] = Encoded_sequences.tolist()
    return dataframe


def pad_values(dataframe, column_name):
    values = dataframe[column_name].to_numpy()
    padded_input = keras.preprocessing.sequence.pad_sequences(
        values,
        maxlen=50,
        padding='post',
        truncating='post',
        dtype='float32'
    )
    
    dataframe[column_name] = padded_input.tolist()
    dataframe[column_name] = dataframe[column_name].apply(lambda x: [(round(value, 2)) for value in x])
    return dataframe


def feature_engineering(dataframe):
    print("Starting with Feature engineering")
    
    dataframe = peptide_length(dataframe)
    print("✓ Completed peptide_length")
    
    dataframe = tryptic(dataframe)
    print("✓ Completed tryptic")
    
    dataframe["Tryptic"] = dataframe["Tryptic"].replace([True, False], [1, 0])
    print("✓ Completed tryptic boolean conversion")
    
    dataframe = basic(dataframe)
    print("✓ Completed basic")
    
    dataframe = acid(dataframe)
    print("✓ Completed acid")
    
    dataframe = non_polar(dataframe)
    print("✓ Completed non_polar")
    
    dataframe = polar(dataframe)
    print("✓ Completed polar")
    
    dataframe = polar_basic(dataframe)
    print("✓ Completed polar_basic")
    
    dataframe = polar_acidic(dataframe)
    print("✓ Completed polar_acidic")
    
    dataframe = aromatic(dataframe)
    print("✓ Completed aromatic")
    
    dataframe = alpha(dataframe)
    print("✓ Completed alpha")
    
    dataframe = beta(dataframe)
    print("✓ Completed beta")
    
    dataframe = turn(dataframe)
    print("✓ Completed turn")
    
    dataframe = missed_cleavages(dataframe)
    print("✓ Completed missed_cleavages")
    
    dataframe = isoelectric_point(dataframe)
    print("✓ Completed isoelectric_point")
    
    dataframe = average_isoelectric_point(dataframe)
    print("✓ Completed average_isoelectric_point")
    
    dataframe = hydrophobicity(dataframe)
    print("✓ Completed hydrophobicity")
    
    dataframe = average_hydrophobicity(dataframe)
    print("✓ Completed average_hydrophobicity")
    
    dataframe = sequence_onehot(dataframe)
    print("✓ Completed sequence_onehot")
    
    dataframe = pad_values(dataframe, "Isoelectric_point_AA")
    print("✓ Completed pad_values for Isoelectric_point_AA")
    
    dataframe = pad_values(dataframe, "Hydrophobicity_AA")
    print("✓ Completed pad_values for Hydrophobicity_AA")
    
    dataframe = N_term_AA(dataframe)
    print("✓ Completed N_term_AA")
    
    dataframe = C_term_AA(dataframe)
    print("✓ Completed C_term_AA")
    
    dataframe = N_term_AA_pI(dataframe)
    print("✓ Completed N_term_AA_pI")
    
    dataframe = N_term_AA_pI_last_2(dataframe)
    print("✓ Completed N_term_AA_pI_last_2")
    
    dataframe = C_term_AA_pI(dataframe)
    print("✓ Completed C_term_AA_pI")
    
    dataframe = C_term_AA_pI_last_2(dataframe)
    print("✓ Completed C_term_AA_pI_last_2")
    
    dataframe = elemental_comp(dataframe)
    print("✓ Completed elemental_comp")
    
    print("Starting elemental_comp_modifications (this may take a while)...")
    dataframe = elemental_comp_modifications(dataframe)
    print("✓ Completed elemental_comp_modifications")
    
    dataframe = monoisotopic_mass(dataframe)
    print("✓ Completed monoisotopic_mass")
    
    dataframe = average_mass(dataframe)
    print("✓ Completed average_mass")
    
    print("Finished with Feature engineering")
    return dataframe
