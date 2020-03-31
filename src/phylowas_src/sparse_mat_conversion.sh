#!/bin/bash
#SBATCH --job-name=get_biomarkers
#SBATCH --partition=dpwall
#SBATCH --mem=100GB
#SBATCH --output=/oak/stanford/groups/dpwall/users/briannac/logs/sparse_mat.out
#SBATCH --error=/oak/stanford/groups/dpwall/users/briannac/logs/sparse_mat.err
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36

python3 /oak/stanford/groups/dpwall/users/briannac/phyloWAS/scripts/sparse_mat_conversion.py

