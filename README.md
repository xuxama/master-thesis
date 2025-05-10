# Peptide selection pipeline

A comprehensive pipeline for selecting peptides for targeted mass spectrometry, integrating multiple prediction tools and filters.

# Overview

This pipeline integrates several tools to facilitate the selection of optimal peptides for targeted mass spectrometry experiments:

- In silico digestion
- Basic property analysis (peptide length and methionine content)
- Uniqueness evaluation with Unipept
- Detectability prediction with DeepMSPeptide and DeepDetect
- Charge state prediction with CPred
- MS2 spectra prediction with MS²PIP

# Environment setup

This pipeline uses multiple conda environments to manage tool dependencies. All environment files are provided in the `environments/` directory.

## Creating environments

```bash
# Create and activate environments from YAML files
conda env create -f environments/digest_env.yml
conda env create -f environments/unipept_env.yml
conda env create -f environments/deepms_env.yml
conda env create -f environments/deepdetect_env.yml
conda env create -f environments/cpred_env.yml
conda env create -f environments/ms2pip_env.yml
```

## Environment descriptions

| Environment | Purpose | Key Dependencies |
| --- | --- | --- |
| digest_env | Basic operations and file conversions | Python, pandas, numpy |
| unipept_env | Unipept peptide uniqueness analysis | Unipept CLI tools |
| deepms_env | DeepMSPeptide predictions | Python 3.6, TensorFlow 1.13.1 |
| deepdetect_env | DeepDetect predictions | Python 3.6, TensorFlow 1.10.0 |
| cpred_env | CPred charge state predictions | Python 3.9, TensorFlow 2.15.1 |
| ms2pip_env | MS²PIP spectra predictions | Python 3.9, ms2pip 4.1.0 |

# Before running

When using the pipeline multiple times, it's recommended to remove all input and output files from previous runs to avoid issues with file overwriting.

# Pipeline steps

## 1. In silico digest

*Skip this step if starting from a list of peptides of interest*

This step performs in silico digestion of proteins to generate peptides.

```bash
# Activate environment
conda activate digest_env

# Execute for a single protein
python in_silico_digest.py <UniProtID> [missed_cleavages]

# OR for multiple proteins in a FASTA file
python in_silico_digest_fasta.py file.fasta [missed_cleavages]
```

**Output:** `peptides.csv`

## 2. Basic properties analysis

This step calculates peptide length and methionine content.

If step 1 was omitted, first create a peptides.csv file containing:

- At minimum, a column named ‘peptide’ with peptide sequences
- Optionally,  a column nmaed ‘Protein_ID’ with corresponding parent sequence identifiers

Example **peptides.csv**

```bash
Protein_ID,peptide
P48632,ETPQSITVVTR
Q9HWL3,LNFNNLFDK
...,...
```

```bash
# Activate environment
conda activate digest_env

# Execute
python basic_properties.py
```

Output: `peptides.csv` with added columns ‘length’ and ‘M_count’

## 3. Unipept peptide uniqueness analysis

*Can be skipped if uniqueness is not a concern*

### 3.1 Create txt file

```bash
# Activate environment
conda activate digest_env

# Execute
python peptides2txt.py peptides.csv

# Deactivate environment
conda deactivate
```

Output: `peptides.txt` containing all unique peptides from peptides.csv

### 3.2 Run Unipept

```bash
# Activate environment
conda activate unipept_env

# Execute Unipept CLI (may take some time)
unipept pept2prot --output peptides_unipept.txt --select peptide,uniprot_id,protein_name,taxon_id --input peptides.txt
```

**Output:** `peptides_unipept.txt` with Unipept results

### 3.2 Process Unipept results

```bash
# Activate environment
conda activate digest_env

# Execute (adjust parameters as needed)
python retrieve_unipept.py --taxon_id 208964 --protein_id P48632 --keywords "Ferripyoverdine receptor,Ferripyoverdine,FpvA" --unipept_file peptides_unipept.txt --peptides_file peptides.csv
```

**Details:**

- Keywords are used to identify matches corresponding to the protein of interest
- The script counts UniProt matches that correspond to the given taxon ID but not to the specified protein_ID or keywords
- Adds columns 'unique' (boolean) and 'problematic_matches' in case of non-unique to `peptides.csv`
- Peptides not found in Unipept results will have missing values

## 4. DeepMSPeptide prediction

### 4.1 Create txt file

*Skip if step 3 was completed*

```bash
# Activate environment
conda activate digest_env

# Execute
python peptides2txt.py peptides.csv

# Deactivate environment
conda deactivate
```

**Output:** `peptides.txt` with unique sequences

### 4.2 Run DeepMSPeptide

```bash
# Copy text file to DeepMSPeptide directory
cp peptides.txt ../deepmspeptide/DeepMSPeptide/DeepMSPeptide

# Activate environment
conda activate deepms_env

# Navigate to DeepMSPeptide directory
cd ../deepmspeptide/DeepMSPeptide/DeepMSPeptide

# Execute
python DeepMSPeptide.py peptides.txt

# Rename results file
mv peptides_Predictions.txt peptides_DMP.txt

# Copy results back to pipeline directory
cp peptides_DMP.txt ../../../peptide_selection_pipeline
```

**Output:** `peptides_DMP.txt` with DeepMSPeptide predictions

### 4.3 Process DeepMSPeptide results

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python retrieve_dmp_preds.py
```

**Output:** Updates `peptides.csv` with added column 'DMP_prob'

## 5. DeepDetect

### 5.1 Run DeepDetect

```bash
# Activate environment
conda activate deepdetect_env

# Navigate to DeepDetect directory
cd deepdetect/SourceCode/deepdetect_pred

# Download FASTA if needed
wget https://www.uniprot.org/uniprotkb/P48632.fasta

# Execute (adjust parameters as needed)
python main.py --input=P48632.fasta --output=peptides_DD.txt --regular='>(\S+)' --protease='Trypsin' --missed_cleavages=2 --min_len=8 --max_len=25

# Copy results to pipeline directory
cp peptides_DD.txt ../../../peptide_selection_pipeline
```

**Output:** `peptides_DD.txt` with DeepDetect results

### 5.2 Process DeepDetect results

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python retrieve_dd_preds.py
```

**Output:** Updates `peptides.csv` with added column 'DD_prob'

## 6. CPred charge state prediction

### 6.1 Prepare CPred input

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python peptides2cpredinput.py
```

**Output:** `cpred_input.csv` formatted for CPred

### 6.2 Run CPred

```bash
# Activate environment
conda activate cpred_optimized_env

# Navigate to CPred directory
cd cpred_optimized/CPred/CPred

# Execute Feature Engineering
python CPred_main.py FeatureEngineering -i ../../../peptide_selection_pipeline/cpred_input.csv -o output_FE

# Execute prediction
python CPred_main.py prediction -i ../output_FE.parquet -m Data/Models/CPred_model_v1.h5 -o predictions

# Rename predictions file
cd predictions
mv Model_predictions.csv peptides_CPred.csv

# Copy to pipeline directory
cp peptides_CPred.csv ../../../../peptide_selection_pipeline
```

**Output:** `peptides_CPred.csv` with CPred results

### 6.3 Process CPred results

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python retrieve_cpred_preds.py
```

**Output:** `peptides.csv` with added columns 'CS' (charge state) and 'CS_prob' (corresponding probability)

## 7. MS²PIP prediction

### 7.1 Prepare MS²PIP input

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python peptides2ms2pipinput.py
```

**Output:** `ms2pip_input.tsv` formatted for MS²PIP

### 7.2 Run MS²PIP

```bash
# Activate environment
conda activate ms2pip_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
ms2pip predict-batch ms2pip_input.tsv -f spectronaut -o peptides_ms2pip --model HCD
```

**Output:** `peptides_ms2pip.spectronaut.tsv` with MS²PIP predictions

### 7.3 Process MS²PIP results

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python retrieve_ms2pip_preds.py
```

**Output:** `peptides.csv` with added columns 'valid_peaks' and 'average_peak_intensity'

## 8. Final peptide selection

```bash
# Activate environment
conda activate digest_env

# Navigate to pipeline directory
cd peptide_selection_pipeline

# Execute
python peptide_selection.py
```

**Output:** `final_selection.csv` containing peptides that meet all filter criteria
