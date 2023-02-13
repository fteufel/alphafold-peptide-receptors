'''
Script to extract ranked lists from checkpoint dir as done in the notebook.

'''
import pandas as pd
import os
from qc_metrics import extract_af_metrics
from concurrent.futures import ProcessPoolExecutor
from tqdm.auto import tqdm
import numpy as np

RECEPTOR_CSV_PATH = 'data/human_receptors.csv'
PREDICTION_DIR = 'predictions_smorfs'
RESULTS_CSV_PATH = 'results/scores_smorfs.csv'


df_receptors = pd.read_csv(RECEPTOR_CSV_PATH)
df_receptors = df_receptors.loc[~df_receptors['is_single_tm']]
df_receptors = df_receptors.loc[df_receptors['Sequence'].str.len()<=2000]
df_receptors.shape

# list all predictions that are complete
predictions = os.listdir(PREDICTION_DIR)
print('Filtering for completed runs only...')
predictions = [x for x in predictions if 'timings.json' in os.listdir(os.path.join(PREDICTION_DIR, x))]

if os.path.exists(RESULTS_CSV_PATH):
    print('Found existing results, filtering already done runs...')
    df =  pd.read_csv(RESULTS_CSV_PATH)
    predictions = [x for x in predictions if x not in df['pred_folder'].values]


n_jobs = 40
with ProcessPoolExecutor(n_jobs) as executor:
    jobs = {}
    for f in tqdm(predictions, leave=False):
        p = os.path.join(PREDICTION_DIR, f)
        file_protein_id = p.split('/')[-1].split('_')[0]
        if file_protein_id not in df_receptors['protein_id'].values:
            continue

        topology = df_receptors.loc[df_receptors['protein_id']==file_protein_id].iloc[0]['Topology']

        job = executor.submit(extract_af_metrics,p, receptor_chains = [0], ligand_chains = [1], receptor_topology_string=topology, relaxed=False, pred_is_nested=False)

        jobs[p] = job

    failed_runs = []
    res = []
    for p, job in tqdm(jobs.items(), leave=True):
        if job.exception() is not None:
            print(p)
            print(job.exception())
            continue
            
        df = job.result()
        if len(df) == 0:
            failed_runs.append(p)
            continue

        s = df.reset_index() # make a new column that contains 1,2,3,4,5 model num
        s['pred_folder'] = df['pred_folder'][0].split('/')[-1]
            
        res.append(s)

df = pd.concat(res)

tmp = df['pred_folder'].str.split('_', expand=True)
df['receptor'] = tmp[0]
df['peptide'] = tmp[1]
df = df.set_index(['peptide', 'receptor', 'index'])

# uupdate previous results if existing
if os.path.exists(RESULTS_CSV_PATH):
    df_old =  pd.read_csv(RESULTS_CSV_PATH)
    df_old = df_old.set_index(['peptide', 'receptor', 'index'])
    df = pd.concat([df_old, df])

df.to_csv(RESULTS_CSV_PATH)


def compute_mad(x):
    return np.median(np.absolute(x - np.median(x)))

df_mad = df.groupby(level=[0,1]).agg(compute_mad)
df_median = df.groupby(level=[0,1]).agg('median')

df_lowerbound = df_median - df_mad
df_lowerbound['PAE'] = df_median['PAE'] + df_mad['PAE'] #PAE low is better
df_lowerbound['contacts_outside'] = df_median['contacts_outside']
df_lowerbound['contacts_inside'] =df_median['contacts_inside']



# enrich output
import requests
annotations = {}
for protein_id in df_lowerbound.index.get_level_values(1).unique().values:
    data = requests.get(f"https://rest.uniprot.org/uniprotkb/{protein_id}.json").json()
    try:
        annotations[protein_id] = [data['proteinDescription']['recommendedName']['fullName']['value'], data['sequence']['length']]
    except KeyError: #will not work for trembl
        pass
#https://rest.uniprot.org/uniprotkb/P03979.json
df_annot = pd.DataFrame.from_dict(annotations).T
df_annot.columns = ['Full name', 'length']


# import ipdb; ipdb.set_trace()
#dumb but works.
df_lowerbound = df_lowerbound.reset_index().set_index('receptor').join(df_annot).reset_index().set_index(['peptide', 'index'])
df_lowerbound.index.names = ['peptide', 'receptor']


excel_dir = RESULTS_CSV_PATH[:-4]
os.makedirs(excel_dir, exist_ok=True)
for pep_id, df_pep in df_lowerbound.groupby(level=0):
    df_pep.sort_values('iptm', ascending=False).to_excel(os.path.join(excel_dir, f'{pep_id.replace(":", "_")}_human_canonical.xlsx'))
