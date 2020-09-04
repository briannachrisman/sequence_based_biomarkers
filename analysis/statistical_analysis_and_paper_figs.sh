#!/bin/bash
#SBATCH --job-name=statistics_asd
#SBATCH --partition=bigmem
#SBATCH --array=1%1
#SBATCH --output=/scratch/groups/dpwall/personal/briannac/logs/statistics_asd_save_each_iter%p_%a.out
#SBATCH --error=/scratch/groups/dpwall/personal/briannac/logs/statistics_asd_save_each_iter%p_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=1000GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu


### file at /scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/analysis/statistical_analysis_and_paper_figs.sh
module load python/3.6.1

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36

cd /scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers


python3 -u analysis/statistical_analysis_and_paper_figs.py 1 #$SLURM_ARRAY_TASK_ID



