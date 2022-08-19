'''
Adapted from original repo.

Use compute_pdockq(pdb_file, pickle_file) instead of command line script.
https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
'''
import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict
import pdb

# parser = argparse.ArgumentParser(description = '''Calculate a predicted DockQ score for a predicted structure.''')
# #Bench4
# parser.add_argument('--pdbfile', nargs=1, type= str, default=sys.stdin, help = 'Path to pdbfile to be scored. Note that this file needs to contain at least two chains')
# parser.add_argument('--pickle_file', nargs=1, type= str, default=sys.stdin, help = 'Path to .pkl file from AlphaFold containing plddt scores.')


#####################FUNCTIONS#########################
def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record

def read_pdb(pdbfile):
    '''Read a pdb file predicted with AF and rewritten to conatin all chains
    '''

    chain_coords = {}
    with open(pdbfile, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            #Get CB - CA for GLY
            if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                if record['chain'] in [*chain_coords.keys()]:
                    chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                else:
                    chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]

    return chain_coords


def calc_pdockq(chain_coords, plddt, t):
    '''Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    '''

    #Get interface
    ch1, ch2 = [*chain_coords.keys()]
    coords1, coords2 = np.array(chain_coords[ch1]), np.array(chain_coords[ch2])
    #Check total length
    if coords1.shape[0]+coords2.shape[0]!=plddt.shape[0]:
        print('Length mismatch with plDDT:',  coords1.shape[0]+coords2.shape[0], plddt.shape[0])

    #Calc 2-norm
    mat = np.append(coords1, coords2,axis=0)
    a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1,l1:]
    contacts = np.argwhere(contact_dists<=t)

    if contacts.shape[0]<1:
        pdockq=0
        ppv=0
    else:
        #Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt[np.unique(contacts[:,0])], plddt[np.unique(contacts[:,1])]]))
        #Get the number of interface contacts
        n_if_contacts = contacts.shape[0]

        x = avg_if_plddt*np.log10(n_if_contacts+1) #Add 1 to avoid NaNs
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

        #PPV
        PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,
            0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
            0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,
            0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
            0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
            0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
            0.63555449, 0.55890174])

        pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
            0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,
            0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,
            0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,
            0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
            0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
            0.06968505, 0.02946438])
        inds = np.argwhere(pdockq_thresholds>=pdockq)
        if len(inds)>0:
            ppv = PPV[inds[-1]][0]
        else:
            ppv = PPV[0]

    return pdockq, ppv


def compute_pdockq(pdb_file, pickle_file):
    #Read chains
    chain_coords = read_pdb(pdb_file)
    #Check chains
    if len(chain_coords.keys())<2:
        print('Only one chain in pdbfile', pdb_file)
        sys.exit()
    #Get plDDT
    pickle_data = np.load(pickle_file, allow_pickle=True)
    plddt = pickle_data['plddt']
    #Calculate pdockq
    t=8 #Distance threshold, set to 8 Å
    pdockq, ppv = calc_pdockq(chain_coords, plddt, t)

    return pdockq



#################MAIN####################

# #Parse args
# args = parser.parse_args()
# #Read chains
# chain_coords = read_pdb(args.pdbfile[0])
# #Check chains
# if len(chain_coords.keys())<2:
#     print('Only one chain in pdbfile', args.pdbfile[0])
#     sys.exit()
# #Get plDDT
# pickle_data = np.load(args.pickle_file[0], allow_pickle=True)
# plddt = pickle_data['plddt']
# #Calculate pdockq
# t=8 #Distance threshold, set to 8 Å
# pdockq, ppv = calc_pdockq(chain_coords, plddt, t)
# print('pDockQ =',np.round(pdockq,3),'for',args.pdbfile[0])
