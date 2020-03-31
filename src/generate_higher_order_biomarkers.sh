#!/bin/bash
#SBATCH --job-name=get_biomarkers
#SBATCH --partition=dpwall
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --nodes=1
#SBATCH --output=/oak/stanford/groups/dpwall/users/briannac/logs/get_biomarkers.out
#SBATCH --error=/oak/stanford/groups/dpwall/users/briannac/logs/get_biomarkers.err
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=briannac@stanford.edu

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36
 
for study in obese_lean_twins diabetes_150 crohns
do
for msa_method in DECIPHER 
do
for order in 1 2
do
echo running $study $msa_method $order

#python3 /oak/stanford/groups/dpwall/users/briannac/phyloWAS/scripts/generate_higher_order_biomarkers.py $order #/oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/${study}/tables/person_vs_taxa.csv #/oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/${study}/tables/${msa_method}_taxa_vs_variants.csv #/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/${study}/person_biomarker/${msa_method}/

done
done
done



for study in yogurt
do
for msa_method in DECIPHER 
do
for order in 3
do
echo running $study $msa_method $order

python3 /oak/stanford/groups/dpwall/users/briannac/phyloWAS/scripts/generate_higher_order_biomarkers.py $order /oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/${study}/tables/person_vs_taxa.csv /oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/${study}/tables/${msa_method}_taxa_vs_variants.csv /oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/${study}/person_biomarker/${msa_method}/

done
done
done


