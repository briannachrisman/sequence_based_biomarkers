#!/bin/bash
#SBATCH --job-name=get_biomarkers
#SBATCH --partition=bigmem
#SBATCH --output=/scratch/groups/dpwall/personal/briannac/logs/get_biomarkers.out
#SBATCH --error=/scratch/groups/dpwall/personal/briannac/logs/get_biomarkers.err
#SBATCH --time=20:00:00
#SBATCH --mem=1000GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu


### file at /scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/src/generate_higher_order_biomarkers.sh
module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36

cd /scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers

for study in yogurt
do
for msa_method in DECIPHER
do
for order in 3
do
echo running $study $msa_method $order
echo "removing duplicates"
#python3 src/remove_duplicate_variants.py $study ${msa_method}
echo "geneating higher order biomarkers"
#python3 src/generate_higher_order_biomarkers.py $order data/${study}/abundance.tsv data/${study}/${msa_method}_taxa_vs_variants_unique.tsv results/person_biomarker/${study}_${msa_method}_
echo "removing SBBs in LD"
python3 -u src/remove_SBB_in_LD.py $study ${msa_method} $order

done
done
done





