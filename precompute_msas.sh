#! /bin/bash
# A name for the job
#SBATCH --job-name=msas
 
# let's try the high memory nodes
#SBATCH --partition=compute
#SBATCH --ntasks=1 --cpus-per-task=12 --mem=10G

#SBATCH --time=48:00:00
#SBATCH --output=/nfs_home/users/fegt/alphafold-peptide-receptors/logs/%j.out
#SBATCH --error=/nfs_home/users/fegt/alphafold-peptide-receptors/logs/%j.err


cd /nfs_home/users/fegt/alphafold-2.2.0
source /nfs_home/users/fegt/.bashrc
source /nfs_home/users/fegt/miniconda3/etc/profile.d/conda.sh

conda activate alphafold
python3 ../alphafold-peptide-receptors/af_scripts/precompute_msas.py