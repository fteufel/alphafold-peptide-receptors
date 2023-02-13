#! /bin/bash
# A name for the job
#SBATCH --job-name=af
#SBATCH --partition=compute
#SBATCH --ntasks=1 --cpus-per-task=32 --mem=50G

#SBATCH --time=128:00:00
#SBATCH --output=/nfs_home/users/fegt/alphafold-peptide-receptors/logs/%j.out
#SBATCH --error=/nfs_home/users/fegt/alphafold-peptide-receptors/logs/%j.err

cd /nfs_home/users/fegt/alphafold-peptide-receptors
source /nfs_home/users/fegt/.bashrc
source /nfs_home/users/fegt/miniconda3/etc/profile.d/conda.sh

conda activate alphafold
python3 gather_results.py

