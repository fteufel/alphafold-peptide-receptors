'''
Start alphafold jobs from precomputed MSAs.
Execute this while in the alphafold-2.2.0 directory.
'''
import subprocess
import argparse
from Bio import SeqIO
import os
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from tqdm.auto import tqdm
from time import sleep
import shutil
import json
import pandas as pd
import pickle

from sqlite_job_manager import SQLiteJobManager

#AF_CONFIG_STR = '--data_dir=../weights --model_preset=multimer --num_multimer_predictions_per_model=1 --max_template_date=1950-11-01 --run_relax=False --uniref90_database_path=../weights/uniref90/uniref90.fasta --mgnify_database_path=../weights/mgnify/mgy_clusters_2018_12.fa --template_mmcif_dir=../weights/pdb_mmcif/mmcif_files --obsolete_pdbs_path=../weights/pdb_mmcif/obsolete.dat --db_preset=reduced_dbs --small_bfd_database_path=../weights/small_bfd/bfd-first_non_consensus_sequences.fasta --pdb_seqres_database_path=../weights/pdb_seqres/pdb_seqres.txt --uniprot_database_path=../weights/uniprot/uniprot.fasta --use_gpu_relax=False --use_precomputed_msas=True'
AF_CONFIG_STR = '--data_dir=../alphafold_data/weights --model_preset=multimer --num_multimer_predictions_per_model=1 --max_template_date=1950-11-01 --run_relax=False --uniref90_database_path=../alphafold_data/weights/uniref90/uniref90.fasta --mgnify_database_path=../alphafold_data/weights/mgnify/mgy_clusters_2022_05.fa --template_mmcif_dir=../alphafold_data/weights/pdb_mmcif/from_einstein/mmcif_files --obsolete_pdbs_path=../alphafold_data/weights/pdb_mmcif/from_einstein/obsolete.dat --db_preset=reduced_dbs --small_bfd_database_path=../alphafold_data/weights/small_bfd/bfd-first_non_consensus_sequences.fasta --pdb_seqres_database_path=../alphafold_data/weights/pdb_mmcif/from_einstein/pdb_seqres.txt --uniprot_database_path=../alphafold_data/weights/uniprot/uniprot.fasta --use_gpu_relax=False --use_precomputed_msas=True'



# Need a lock for thread-safe GPU availability management
# LOCK = Lock()
# # We use this dict for the ThreadPool jobs to block a GPU, so that a new job will not launch on the same one.
# # NOTE adjust number of gpus/gpus reserved for use manually
# GPU_AVAILABLE = {
#     0 : True,
#     # 1 : True,
#     # 2 : True,
#     # 3 : True,
#     # 4 : True,
#     # 5 : True,
#     # 6 : True,
#     # 7 : True,
# }

# def get_gpu():
#     '''Get a GPU and block it. If none available, wait and retry.'''
#     while True:
#         with LOCK:
#             for k, v in GPU_AVAILABLE.items():
#                 if v:
#                     GPU_AVAILABLE[k] = False
#                     return k

#         print('No GPU available, retry later.')
#         sleep(60)

# def unblock_gpu(gpu_id):
#     with LOCK:
#         GPU_AVAILABLE[gpu_id] = True
    
def cleanup_pickle(fp):
    '''Remove some data from the pickles to reduce disk usage.'''
    def try_delete_key(d, k):
        '''Don't fail when trying to re-delete.'''
        try:
            del d[k]
        except KeyError:
            pass

    # reduce the pickle size
    output = pickle.load(open(fp, 'rb'))

    try_delete_key(output, 'masked_msa')
    try_delete_key(output, 'experimentally_resolved')
    try_delete_key(output, 'structure_module')
    try_delete_key(output, 'predicted_lddt')
    try_delete_key(output, 'distogram')
    try_delete_key(output, 'aligned_confidence_probs')

    pickle.dump(output, open(fp, 'wb'))


def assemble_msa(msa_dir, chain_1_name, chain_2_name, chain_1_seq, chain_2_seq, output_dir):
    '''Assembles two checkpointed MSAs into the multimer input format.'''

    os.makedirs(output_dir, exist_ok=True)

    msa_path = os.path.join(msa_dir, chain_1_name)
    shutil.copytree(msa_path, os.path.join(output_dir ,'msas', 'A'), dirs_exist_ok=True)
    msa_path = os.path.join(msa_dir, chain_2_name)
    shutil.copytree(msa_path, os.path.join(output_dir, 'msas', 'B'), dirs_exist_ok=True)

    # make the id mapping json
    chain_id_map = {
    "A": {
        "description": chain_1_name,
        "sequence": chain_1_seq,
    },
    "B": {
        "description": chain_2_name,
        "sequence": chain_2_seq,
    }
    }
    json.dump(chain_id_map, open(os.path.join(output_dir, 'msas', 'chain_id_map.json'), 'w'), indent=4)

    return 



def alphafold_job(manager, msa_dir, log_dir):

    pep_id, rec_id, pep_seq, rec_seq, out_dir = manager.get_job()

    # we just do nothing if the database is empty.
    if pep_id is None:
        return None


    # rec_id, pep_id, rec_seq, pep_seq, msa_dir, log_dir = x
    complex_seqs = f'{rec_id}_{pep_id}.fasta'
    open(complex_seqs, 'w').write(
f""">{rec_id}
{rec_seq}
>{pep_id}
{pep_seq}
"""
                )

    # gpu_id = get_gpu()

    try:
        assemble_msa(msa_dir, rec_id, pep_id, rec_seq, pep_seq, os.path.join(out_dir, f'{rec_id}_{pep_id}'))

        my_env = os.environ.copy()
        # my_env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

        af_command = f'python3 run_alphafold.py --fasta_paths={complex_seqs} --output_dir {out_dir} {AF_CONFIG_STR}'
        with open(os.path.join(log_dir, f'{complex_seqs}.log'), 'w') as f:
            subprocess.run(af_command, env=my_env, shell=True, stdout=f, stderr=subprocess.STDOUT) #redirect output to a log file for debugging.
        # unblock_gpu(gpu_id)
        
        # clean up inputs and outputs
        shutil.rmtree(os.path.join(out_dir, f'{rec_id}_{pep_id}', 'msas'))
        os.remove(os.path.join(out_dir, f'{rec_id}_{pep_id}', 'features.pkl'))

        pickles = [x for x in os.listdir(os.path.join(out_dir, f'{rec_id}_{pep_id}')) if 'result_model_' in x]
        for p in pickles:
            cleanup_pickle(os.path.join(out_dir, f'{rec_id}_{pep_id}', p))

        os.remove(complex_seqs)

        manager.mark_job_as_complete(pep_id, rec_id, out_dir)
    except Exception as e:
        print(e)
        # unblock_gpu(gpu_id)
        manager.mark_job_as_failed(pep_id, rec_id, out_dir)

        
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--peptides', required=True, help='Fasta file of peptides to run.')
    parser.add_argument('--receptors', required=True, help='DeepTMHMM annotated receptor library in .csv format.')
    parser.add_argument('--out_dir', required=True, help='Output directory for alphafold predictions.')
    parser.add_argument('--msa_dir', required=True, help='Directory with precomputed MSAs.')
    parser.add_argument('--min_seq_len', default=1)
    parser.add_argument('--max_seq_len', default=2000)
    parser.add_argument('--log_dir', default='logs')
    parser.add_argument('--max_jobs', type=int, default=100)
    parser.add_argument('--job_db', type=str, default = 'alphafold_jobs.sqlite')
    parser.add_argument('--skip_db_update', action='store_true', help='Do not process input files + out dir and try adding missing to db.')

    # PEPTIDE_FASTA = '../data/pdb_benchmark_peptides.fasta'
    # RECEPTOR_CSV = '../data/human_receptors.csv'
    # OUT_DIR = '../predictions/'
    # MSA_DIR = '../data/msas'
    # MAX_SEQ_LEN = 2000
    # MIN_SEQ_LEN = 1
    # NUM_PARALLEL = 8
    # LOG_DIR = '../logs'
    args = parser.parse_args()


    # filter receptor list and create missing MSAs.
    df = pd.read_csv(args.receptors)#, index_col=[0,1])

    
    df = df.loc[df['Sequence'].str.len()<=args.max_seq_len]
    df = df.loc[df['Sequence'].str.len()>=args.min_seq_len]




           
    # predict structures
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)


    manager = SQLiteJobManager(args.job_db)
    # manager.reset_failed_jobs()
    
    if not args.skip_db_update:
    # inputs = []
        already_predicted = os.listdir(args.out_dir)
        already_predicted = [x for x in already_predicted if 'timings.json' in os.listdir(os.path.join(args.out_dir, x))]

        for pep in SeqIO.parse(open(args.peptides), 'fasta'):
            _, prot_id, pos = pep.name.split('|')
            pep_id = prot_id + ':' + pos
            pep_seq = str(pep.seq)


            for idx, row in df.iterrows():
                rec_id = row['protein_id']
                rec_seq = row['Sequence']
                
                if len(rec_seq)>args.max_seq_len:
                    print(f'Skipping {rec_id}: {len(rec_seq)} AAs.')
                    continue
                elif f'{rec_id}_{pep_id}' in already_predicted:
                    print(f'Already processed {rec_id}_{pep_id} .')
                    continue
                else:
                    manager.try_add_job(pep_id, rec_id, pep_seq, rec_seq, args.out_dir)
                    # inputs.append((rec_id, pep_id, rec_seq, pep_seq))

    for i in tqdm(range(args.max_jobs), total=args.max_jobs):
        r = alphafold_job(manager, args.msa_dir, args.log_dir)

    # with ThreadPoolExecutor(4) as executor:
    #     jobs = []
    #     # for inp in inputs:
    #     for i in range(args.max_jobs):
    #         j = executor.submit(alphafold_job, manager, args.msa_dir)
    #         jobs.append(j)

    #     for j in tqdm(jobs, total=len(jobs)):
    #         j.result()






if __name__ == '__main__':
    main()