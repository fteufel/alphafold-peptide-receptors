import mdtraj as md
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
import glob
from pdockq import compute_pdockq



def extract_af_metrics(pred_folder, 
                       receptor_chains,
                       ligand_chains,
                       receptor_topology_string = None,
                       relaxed=True,
                       pred_is_nested=False):
    
    df = pd.DataFrame({})
    relaxed_str = 'relaxed' if relaxed else 'unrelaxed'

    # depending on the dir structure, we need different globbing
    # True  pred_folder/XX/*_model_*.pdb
    # False pred_folder/*_model_*.pdb
    if pred_is_nested:
        pdb_files = [f for f in glob.glob(os.path.join(pred_folder, f'*/{relaxed_str}_model_*.pdb'))]
    else:
        pdb_files = [f for f in glob.glob(os.path.join(pred_folder, f'*{relaxed_str}_model_*.pdb'))]
    pdb_files = np.asarray(pdb_files)
    

    for pdb_file in pdb_files:
        model_num = pdb_file.split('/')[-1].split('_')[2]
        df.at[model_num, 'pred_folder'] = pred_folder
        
        # Open pickle file
        pickle_file = [f for f in glob.glob(os.path.join( os.path.dirname(pdb_file), 'result_model_' + str(model_num) + '*.pkl'))][0]
        prediction = pd.read_pickle(pickle_file)

        df.at[model_num, 'pdockq'] = compute_pdockq(pdb_file, pickle_file)
        
        # Extract ptm, iptm and ranking confidence
        df.at[model_num, 'ptm'] = prediction['ptm']
        df.at[model_num, 'iptm'] = prediction['iptm']
        df.at[model_num, 'ranking_confidence'] = prediction['ranking_confidence']
        
        # Extract plddt and PAE average over binding interface
        model_mdtraj = md.load(pdb_file)
        table, bonds = model_mdtraj.topology.to_dataframe()
        table = table[(table['name']=='CA')]
        table['residue'] = np.arange(0, len(table))
        receptor_res = table[table['chainID'].isin(receptor_chains)]['residue']
        ligand_res = table[table['chainID'].isin(ligand_chains)]['residue']
        
        input_to_calc_contacts = []
        for i in ligand_res:
            for j in receptor_res:
                input_to_calc_contacts.append([i,j])
        
        contacts, input_to_calc_contacts = md.compute_contacts(model_mdtraj, contacts=input_to_calc_contacts, scheme='closest', periodic=False)
        receptor_res_in_contact = []
        ligand_res_in_contact = []
        
        for i in input_to_calc_contacts[np.where(contacts[0]<0.35)]: # threshold in nm
            ligand_res_in_contact.append(i[0])
            receptor_res_in_contact.append(i[1])
        receptor_res_in_contact, receptor_res_counts = np.unique(np.asarray(receptor_res_in_contact), return_counts=True)
        ligand_res_in_contact, ligand_res_counts = np.unique(np.asarray(ligand_res_in_contact), return_counts=True)
        
        if len(ligand_res_in_contact) > 0:
            df.at[model_num, 'plddt_ligand'] = np.median(prediction['plddt'][ligand_res_in_contact])
            df.at[model_num, 'plddt_receptor'] = np.median(prediction['plddt'][receptor_res_in_contact])
            
            df.at[model_num, 'PAE'] = np.median(prediction['predicted_aligned_error'][receptor_res_in_contact,:][:,ligand_res_in_contact])


            if receptor_topology_string is not None:
                labels = np.array(list(receptor_topology_string)) # make indexable.
                intracellular_contact_count = (labels[receptor_res_in_contact] == 'I').sum()
                extracellular_contact_count = (labels[receptor_res_in_contact] == 'O').sum()
                transmembrane_contact_count = (labels[receptor_res_in_contact] == 'M').sum()
                df.at[model_num, 'contacts_outside'] = extracellular_contact_count
                df.at[model_num, 'contacts_inside'] =  intracellular_contact_count

        
            
    return df

