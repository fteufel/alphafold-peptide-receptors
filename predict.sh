#! /bin/bash
# A name for the job
#SBATCH --job-name=af
 
# let's try the high memory nodes
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1 --cpus-per-task=4 --mem=50G

#SBATCH --time=128:00:00
#SBATCH --output=/nfs_home/users/fegt/alphafold-peptide-receptors/logs/%j.out
#SBATCH --error=/nfs_home/users/fegt/alphafold-peptide-receptors/logs/%j.err

nvidia-smi

cd /nfs_home/users/fegt/alphafold-2.2.0
source /nfs_home/users/fegt/.bashrc
source /nfs_home/users/fegt/miniconda3/etc/profile.d/conda.sh

conda activate alphafold
python3 ../alphafold-peptide-receptors/af_scripts/predict_from_precomputed_multinode.py \
--peptides /nfs_home/users/fegt/alphafold-peptide-receptors/data/smorfs.fasta \
--receptors /nfs_home/users/fegt/alphafold-peptide-receptors/data/human_receptors.csv \
--out_dir /nfs_home/users/fegt/alphafold-peptide-receptors/predictions_smorfs \
--msa_dir /nfs_home/users/fegt/alphafold_data/msas \
--max_jobs 100 \
--log_dir /nfs_home/users/fegt/alphafold-peptide-receptors/logs/ \
--job_db /nfs_home/users/fegt/alphafold-peptide-receptors/smorfjobs.sqlite

