'''
Start alphafold jobs from precomputed MSAs.
Execute this while in the alphafold-2.2.0 directory.
Make sure that run_alphafold_msaonly.py is there.
'''
import os
import shutil
import pandas as pd
from tqdm.auto import tqdm
from multiprocessing import Pool
import subprocess
from Bio import SeqIO

NUM_PARALLEL = 2
RECEPTOR_CSV = '../alphafold-peptide-receptors/data/human_receptors.csv'
PEPTIDE_FASTA = '../alphafold-peptide-receptors/data/smorfs.fasta'
MSA_DIR = '../alphafold_data/msas'
MAX_SEQ_LEN = 2000
AF_CONFIG_STR = '--data_dir=../alphafold_data/weights --model_preset=multimer --num_multimer_predictions_per_model=1 --max_template_date=1950-11-01 --run_relax=False --uniref90_database_path=../alphafold_data/weights/uniref90/uniref90.fasta --mgnify_database_path=../alphafold_data/weights/mgnify/mgy_clusters_2022_05.fa --template_mmcif_dir=../alphafold_data/weights/pdb_mmcif/from_einstein/mmcif_files --obsolete_pdbs_path=../alphafold_data/weights/pdb_mmcif/from_einstein/obsolete.dat --db_preset=reduced_dbs --small_bfd_database_path=../alphafold_data/weights/small_bfd/bfd-first_non_consensus_sequences.fasta --pdb_seqres_database_path=../alphafold_data/weights/pdb_mmcif/from_einstein/pdb_seqres.txt --uniprot_database_path=../alphafold_data/weights/uniprot/uniprot.fasta --use_gpu_relax=False --use_precomputed_msas=True'
OUT_DIR = 'msa_tempdir'
LOG_DIR = '../alphafold-peptide-receptors/logs'


def msa_job(x):
    '''Submit single-seq MSAs also as multimer to make MSA types match.'''
    rec_id, rec_seq = x
    complex_seqs = f'{rec_id}.fasta'
    open(complex_seqs, 'w').write(
f""">{rec_id}
{rec_seq}
>dummypep
AAA
"""
    )

    my_env = os.environ.copy()
    my_env["TMPDIR"] = OUT_DIR

    os.makedirs('../logs', exist_ok=True)
    af_command = f'python3 run_alphafold_msaonly.py --fasta_paths={complex_seqs} --output_dir {OUT_DIR} {AF_CONFIG_STR}'
    with open(os.path.join(LOG_DIR, f'{complex_seqs}.log'), 'w') as f:
        subprocess.run(af_command, env=my_env, shell=True, stdout=f, stderr=subprocess.STDOUT) #redirect output to a log file for debugging.

    
    #
    msa_path = os.path.join(OUT_DIR, rec_id, 'msas', 'A')

    try:
        shutil.copytree(msa_path, os.path.join(MSA_DIR, rec_id), dirs_exist_ok=True)
        shutil.rmtree(os.path.join(OUT_DIR, rec_id))
        os.remove(complex_seqs)
    except FileNotFoundError as e:
        print(e)




if __name__ == '__main__':
    os.makedirs(MSA_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    # filter receptor list and create missing MSAs.
    df = pd.read_csv(RECEPTOR_CSV)
    df = df.loc[~df['is_single_tm']]

    df = df.loc[df['Sequence'].str.len()<=MAX_SEQ_LEN]


    already_predicted = os.listdir(MSA_DIR)
    
    jobs = []
    for idx, row in df.iterrows():

        if row['protein_id'] != 'P35456':
            continue
        seq_name = row['protein_id']
        if seq_name in already_predicted:
            continue
        aa_sequence = row['Sequence']
        jobs.append((seq_name, aa_sequence))

    for pep in SeqIO.parse(open(PEPTIDE_FASTA), 'fasta'):
        _, prot_id, pos = pep.name.split('|')
        pep_id = prot_id + ':' + pos
        pep_seq = str(pep.seq)
        jobs.append((pep_id, pep_seq))


    with Pool(processes=NUM_PARALLEL) as pool:
        for i in tqdm(pool.imap_unordered(msa_job, jobs), total=len(jobs)):
            pass

