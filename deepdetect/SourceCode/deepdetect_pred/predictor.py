from keras.models import model_from_json
from keras.preprocessing.sequence import pad_sequences as padding
import numpy as np
import time
import re
from math import sqrt
import os

# peptide / 31-mer coding
def coding(seq_type, seqs):
    # amino acid dictionaries
    if seq_type == 'peptides':
        dic = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9,
               'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
               'V': 18, 'W': 19, 'Y': 20}
    elif seq_type == '31mers':
        dic = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8,
               'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16,
               'V': 17, 'W': 18, 'Y': 19, 'Z': 20}
    
    # coding
    coded_seqs = []
    for seq in seqs:
        coded_seq = [dic.get(aa) for aa in seq]
        coded_seqs.append(coded_seq)  # Using append instead of += for clarity
    
    return coded_seqs


# predicting
def predictor(protease, data, padding_len, res_path):
    print("Loading model.")
    s3 = time.time()
    
    # loading BiLSTM model with Linux path syntax
    try:
    	print("Attempting to load BiLSTM model from JSON and weights...")
    	json_bilstm = os.path.join('..', '..', 'DeepDetect', f'BiLSTM_{protease}.json')
    	h5_bilstm = os.path.join('..', '..', 'DeepDetect', f'BiLSTM_{protease}.h5')

    	import json
    	import h5py
    	import numpy as np

    	with open(json_bilstm, 'r') as bilstm:
        	model_bilstm = bilstm.read()
    	loaded_model_bilstm = model_from_json(model_bilstm)
    	#Custom weights loading that handles encoding issue
    	if os.path.exists(h5_bilstm):
    		print(f"Loading weights from {h5_bilstm}")
    		try:
    			# Try loading as a full model first
    			from keras.models import load_model
    			temp_model = load_model(h5_bilstm, compile=False)
    			loaded_model_bilstm.set_weights(temp_model.get_weights())
    			print("Successfully loaded weights from full model.")
    		except Exception as e1:
    			print(f"Could not load as full model: {e1}")
    			try:
    				# Try direct weight loading with workaround
    				f = h5py.File(h5_bilstm, 'r')
    				weight_names = [n.decode('utf8') if isinstance(n, bytes) else n
                                	for n in f.attrs.get('weight_names')]
    				if len(weight_names) > 0:
    					weights = []
    					for name in weight_names:
    						weight = np.array(f[name])
    						weights.append(weight)
    					loaded_model_bilstm.set_weights(weights)
    					print("Successfully loaded weights manually.")
    				else:
    					print("No weights found in the H5 file.")
    				f.close()
    			except Exception as e2:
    				print(f"Error with manual weights loading: {e2}")
    				# Fall back to old method as last resort
    				loaded_model_bilstm.load_weights(h5_bilstm)
    	else:
    		print(f"Warning: Weight file {h5_bilstm} not found!")
    	print("BiLSTM model loaded successfully.")
    except Exception as e:
    	print(f"Error loading BiLSTM model: {e}")
    	raise

    # loading DeepDigest model
    json_deepdigest = os.path.join('..', '..', 'DeepDetect', f'DeepDigest_{protease}.json')
    h5_deepdigest = os.path.join('..', '..', 'DeepDetect', f'DeepDigest_{protease}.h5')
    
    with open(json_deepdigest, 'r') as deepdigest:
        model_deepdigest = deepdigest.read()
    loaded_model_deepdigest = model_from_json(model_deepdigest)
    loaded_model_deepdigest.load_weights(h5_deepdigest)
    
    e3 = time.time()
    print("Time cost of loading model is %s seconds." % (e3 - s3))
    
    # prediction
    print("Predicting! Please wait...")
    s4 = time.time()
    
    # extract all the peptides and 31-mers
    peps = []
    mers = []
    for line in data:
        seqs = re.split('[\t,]', line[1])
        peps.append(seqs[0])  # Using append for clarity
        mers.extend(mer for mer in seqs[1:] if len(mer) > 1)  # Using extend instead of += for clarity
    print("There are %s in silico digested peptides.\n"
          "There are %s candidate cleavage sites." % (len(peps), len(mers)))

    #==============================================================================    
    #     BiLSTM prediction
    #==============================================================================    
    coded_peps = coding('peptides', peps)
    x_peps = padding(coded_peps, maxlen=padding_len,
                     padding='post', truncating='post', value=0)
    bilstm_pred = loaded_model_bilstm.predict(x_peps).squeeze()
    
    #==============================================================================    
    #     DeepDigest prediction
    #==============================================================================    
    x_mers = np.array(coding('31mers', mers))
    mers_pred = loaded_model_deepdigest.predict(x_mers).squeeze()
    
    e4 = time.time()
    print("Time cost of prediction is %s seconds." % (e4 - s4))
    
    # =============================================================================
    #     calculating the peptide detectabilities
    # =============================================================================
    print("Calculating peptide detectabilities.")
    s5 = time.time()
    bilstm_ind = 0  # bilstm detectability
    dig_ind = 0  # digestibility
    results = []
    
    for line in data:
        pro_id, seqs = line
        seqs_list = seqs.split('\t')
        bilstm_prob = float(bilstm_pred[bilstm_ind])
        bilstm_ind += 1
        
        # The model predicted the probabilities of missed digestion!!!
        if len(seqs_list) == 4:
            pep, left_mer, right_mer, missed_mers = seqs_list
            
            # digestibility of the left 31-mer
            if left_mer != '*':
                left_dig = 1 - float(mers_pred[dig_ind])
                dig_ind += 1
            else:
                left_dig = float(1)
            
            # digestibility of the right 31-mer
            if right_mer != '*':
                right_dig = 1 - float(mers_pred[dig_ind])
                dig_ind += 1
            else:
                right_dig = float(1)
            
            # digestibilities of the missed 31-mers
            missed_sites = missed_mers.split(',')
            missed_probs = ''
            missed_dig = float(1)
            for site in missed_sites:
                if len(site) == 31:
                    prob = float(mers_pred[dig_ind])
                    dig_ind += 1
                    missed_probs += str(1 - prob) + ','
                    missed_dig *= prob
                elif site != '':
                    print(line)
                    print(site)
                    
            # peptide digestibility
            dig_prob = left_dig * right_dig * missed_dig
            
            # peptide detectability
            det_prob = sqrt(bilstm_prob * dig_prob)
            results.append([pro_id, pep, det_prob])  # Using append for clarity
        elif len(seqs_list) == 1:
            results.append(line + [bilstm_prob])  # Using append for clarity
    
    e5 = time.time()
    print("Time cost of calculation is %s seconds." % (e5 - s5))
    
    # save results
    s6 = time.time()
    np.savetxt(res_path, results, fmt='%s\t%s\t%s',
               delimiter='\t', newline='\n',
               header='Protein id\tPeptide sequence\tPeptide detectability',
               comments='')
    e6 = time.time()
    print("Time cost of saving results is %s seconds." % (e6 - s6))
