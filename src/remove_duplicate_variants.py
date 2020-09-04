### file at /scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/src/remove_duplicate_variants.py

import pandas as pd
import numpy as np
import sys
import csv 

study = sys.argv[1]
msa_method = sys.argv[2]

taxa_vs_variants = pd.read_csv('data/' + study + '/' + str(msa_method) + '_taxa_vs_variants.tsv', sep='\t')
taxa_vs_variants = taxa_vs_variants.transpose().drop_duplicates().transpose()
taxa_vs_variants.columns = ['\"' + c + '\"' for c in taxa_vs_variants.columns]
taxa_vs_variants.to_csv('data/' + study + '/' + str(msa_method) + '_taxa_vs_variants_unique.tsv', sep='\t', quoting=csv.QUOTE_NONE)