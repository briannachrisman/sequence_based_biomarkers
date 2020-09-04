### file at /scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/src/remove_SBB_in_LD.py

import pandas as pd
import numpy as np
import sys

study = sys.argv[1]
msa_method = sys.argv[2]
order = int(sys.argv[3])

print(order, study, msa_method)
biomarkers = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/results/person_biomarker/' + study  + 
                         '_'  + str(msa_method) +  '_biomarkers' + str(order) + '.txt', sep='\t', header=None)
biomarkers = biomarkers.drop(list(range(order)), axis =1)
biomarkers.columns = [i for i in range(order)]
taxa_vs_variants = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/' + str(msa_method) + '_taxa_vs_variants_unique.tsv', sep='\t')
taxa_vs_variants.columns = ['ASV'] + [i for i in range(len(taxa_vs_variants.columns)-1)]
#biomarker_seqs = set()
#biomarkers_seqs_list = []
#unique_idx = [] 

#for i,b in enumerate(biomarkers.iterrows()):
#    if (i%100000)==0:
#        print(i/len(biomarkers), " way done...")
#    biomarker_seq = tuple(sum([taxa_vs_variants[b[1][order_i]].values for order_i in range(order)])==order)
#    if biomarker_seq not in biomarker_seqs:
#        biomarker_seqs = biomarker_seqs.union(set([biomarker_seq]))
#        biomarkers_seqs_list.append(biomarker_seq)
#        unique_idx = unique_idx + [i]
#np.save('results/person_biomarker/' + study  + 
#        '_'  + str(msa_method) +  '_biomarker_seq' + str(order) + '.npy', biomarkers_seqs_list)
#print("step done... deleting biomarker_seq")
#del biomarker_seqs
#del biomarkers_seqs_list

person_variant = np.array(np.load('results/person_biomarker/' + study  + 
                                  '_'  + str(msa_method) +  '_person_variant' + str(order) + '_condensed.npy'))
person_variant, unique_idx = np.unique(person_variant, axis=1, return_index=True)
print(len(unique_idx), " unique markers")
print('here2')
np.save('results/person_biomarker/' + study  + 
        '_'  + str(msa_method) +  '_person_variant' + str(order) + '_unique.npy', person_variant)
del person_variant

biomarkers = biomarkers.iloc[unique_idx]
biomarkers.to_csv('results/person_biomarker/' + study  + 
                                               '_'  + str(msa_method) +  '_biomarkers' + str(order) + '_unique.npy', sep='\t')
print("step done... deleting biomarkers")
del biomarkers
print('here')
person_variant = np.array(np.load('results/person_biomarker/' + study  + 
                                  '_'  + str(msa_method) +  '_person_variant' + str(order) + '_condensed.npy'))
print('here1')
print(len(unique_idx), " unique combos")
print('done, success!')