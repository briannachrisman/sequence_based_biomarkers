import matplotlib
matplotlib.use('Agg')

###############################################
##################### ASD #####################
###############################################

# START
import pandas as pd
import numpy as np
import seaborn as sns
from collections import Counter
import scipy.stats as ss
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from scipy.stats import wilcoxon
from scipy.stats import mannwhitneyu
import multiprocessing
from joblib import Parallel, delayed
import time
import sys

seed = int(time.time())
num_cores = multiprocessing.cpu_count()
print("num cores: ", num_cores)

PROJECT_DIR = '/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers'
PLOT_DIR =  '/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/figs/'
np.random.seed(seed)

FAKE_P_VAL = 2 # This needs to be >1.
N_ITER = 68

study = 'yogurt'
msa_method = 'DECIPHER'
otus = pd.read_csv(PROJECT_DIR + '/data/' + study + '/otu.tsv', sep='\t')
otus.apply(lambda x: x/sum(x), axis=1).to_csv(
    PROJECT_DIR + '/data/' + study + '/abundance.tsv', sep='\t')
tax_table = pd.read_csv(PROJECT_DIR + '/data/' + study + '/tax_table.tsv', sep='\t')

abundance = pd.read_csv(PROJECT_DIR + '/data/' + study + '/abundance.tsv', sep='\t')
person_vs_pheno = pd.read_csv(PROJECT_DIR + '/data/' + study + '/sample_data.tsv', sep='\t')
abundance.index = abundance['Unnamed: 0']
abundance = abundance.drop('Unnamed: 0', axis=1)

person_vs_pheno = pd.read_csv(PROJECT_DIR + '/data/' + study + '/sample_data.tsv', sep='\t')
pairs_counts = Counter([(a,b) for a,b in zip(person_vs_pheno['family'], person_vs_pheno['timepoint'])])
individual_counts = Counter([(a,b,c) for a,b,c in zip(person_vs_pheno['family'], person_vs_pheno['timepoint'], person_vs_pheno['phenotype'])])
good_pairs = [(fam, time) for fam,time in pairs_counts.keys() if (pairs_counts[(fam,time)]==2) & (individual_counts[(fam,time,'A')]==1) & (individual_counts[(fam,time,'N')]==1)]
good_pairs_ids = [(i,j) in good_pairs for i,j in zip(person_vs_pheno['family'], person_vs_pheno['timepoint'])]
person_vs_pheno = person_vs_pheno[good_pairs_ids]
person_vs_pheno = person_vs_pheno.sort_values(['family', 'timepoint', 'phenotype'])
person_vs_pheno['pheno'] = person_vs_pheno['phenotype']=='A'

print("Loading taxa aggregated by category")

# Taxa aggregated by taxa category
person_variant_species = np.array(abundance.transpose().join(tax_table).groupby('Species').aggregate(sum).transpose())
person_variant_genus = np.array(abundance.transpose().join(tax_table).groupby('Genus').aggregate(sum).transpose())
person_variant_family = np.array(abundance.transpose().join(tax_table).groupby('Family').aggregate(sum).transpose())
person_variant_order = np.array(abundance.transpose().join(tax_table).groupby('Order').aggregate(sum).transpose())
person_variant_class = np.array(abundance.transpose().join(tax_table).groupby('Class').aggregate(sum).transpose())
person_variant_phylum = np.array(abundance.transpose().join(tax_table).groupby('Phylum').aggregate(sum).transpose())
person_variant_kingdom = np.array(abundance.transpose().join(tax_table).groupby('Kingdom').aggregate(sum).transpose())


print("Loading taxa aggregated by SBB")
person_variant_1 = np.array(np.load(PROJECT_DIR + '/results/person_biomarker/' + study  + 
                                   '_'  + str(msa_method) +  '_person_variant1_unique.npy'))
person_variant_2 = np.array(np.load(PROJECT_DIR + '/results/person_biomarker/' + study  + 
                            '_'  + str(msa_method) +  '_person_variant2_unique.npy'))
person_variant_3 = np.array(np.load(PROJECT_DIR + '/results/person_biomarker/' + study  + 
                                   '_'  + str(msa_method) +  '_person_variant3_unique.npy'))

tax_table = tax_table.loc[pd.read_csv(PROJECT_DIR + '/data/yogurt/DECIPHER_taxa_vs_variants.tsv', '\t').index] # Fix tax table indexing issue for Yogurt data.


true_phenos = person_vs_pheno['pheno'].values




f_FDR,ax_FDR = plt.subplots(figsize=(10,5))
f_qq,ax_qq = plt.subplots(figsize=(10,5))
f_FDR2,ax_FDR2 = plt.subplots(figsize=(10,5))
for person_variant, title in [#(person_variant_species, 'species'), (person_variant_genus, 'genus'),  (person_variant_family, 'family'),  (person_variant_order, 'order'), (person_variant_class, 'class'), 
                              #(person_variant_1, '1-SBB'), (person_variant_2, '2-SBB'), 
                              (person_variant_3, '3-SBB')]:
    print(title)
    # Compute Statistics
    k = person_vs_pheno['family'].values
    v = [np.random.rand()>.5 for i in person_vs_pheno['family'].values]
    person_variant = person_variant[good_pairs_ids]
    print(np.shape(person_variant)[1], " groups to search through!")
    
    def computeWilcoxon(i, phenos):
        #if (i%10000==0): print(i)
        diffs = person_variant[phenos,i] - person_variant[~phenos,i]
        diffs = diffs[diffs!=0]
        _,p_data = wilcoxon(diffs)
        if np.isnan(p_data):
            return FAKE_P_VAL
        else: 
            return p_data
    
    # Statistics
    #if seed==1:
    if False:
        R_data = np.array(Parallel(n_jobs=num_cores-2)(delayed(computeWilcoxon)(i, true_phenos)
                                                        for i in range(np.shape(person_variant)[1])))
        #print(R_data)
        R_data_na = R_data[R_data!=FAKE_P_VAL]
        np.save(PROJECT_DIR + '/results/sig_values/yogurt_' + 'R_data_' + title + '.npy',R_data_na)

    R_permute = np.zeros((np.shape(person_variant)[1], N_ITER))
    for r in range(47,N_ITER):
        t = time.time()
        print(r, ' iteration')
        v = [np.random.rand()>.5 for i in person_vs_pheno['family'].values]
        to_flip_dict = {i:j for i,j in zip(k,v)}
        new_phenos = np.array(
            [not pheno if to_flip_dict[fam] else pheno for fam, pheno in zip(person_vs_pheno['family'], person_vs_pheno['pheno'])])
        R_permute[:,r] = Parallel(n_jobs=num_cores-2)(delayed(computeWilcoxon)(i, new_phenos)
                                                        for i in range(np.shape(person_variant)[1]))
        print(time.time()-t, ' time taken')
        R_permute_small = R_permute[:,r]
        R_permute_small = R_permute_small[R_permute_small!=FAKE_P_VAL]
        np.save(PROJECT_DIR + '/results/sig_values/yogurt_' + 'R_permute_small' + title + '_iter' + str(r) + '.npy',R_permute_small)
    R_permute = R_permute.flatten()
    R_permute_na = R_permute[R_permute!=FAKE_P_VAL]
    np.save(PROJECT_DIR + '/results/sig_values/yogurt_' + 'R_permute_' + title + str(seed) + '.npy', R_permute_na)

print("success!")
exit()


###############################################
################### Obesity ###################
###############################################

### START HERE ####
study = 'obese_lean_twins'
msa_method = 'DECIPHER'
otus = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/otu.tsv', sep='\t')
otus.apply(lambda x: x/sum(x), axis=1).to_csv(
    '/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/abundance.tsv', sep='\t')
tax_table = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/tax_table.tsv', sep='\t')

abundance = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/abundance.tsv', sep='\t')
person_vs_pheno = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/sample_data.tsv', sep='\t')
abundance.index = abundance['Unnamed: 0']
abundance = abundance.drop('Unnamed: 0', axis=1)

person_vs_pheno = pd.read_csv('/scratch/groups/dpwall/personal/briannac/sequence_based_biomarkers/data/' + study + '/sample_data.tsv', sep='\t')
person_vs_pheno['pheno'] = person_vs_pheno['obesitycat']=='Obese'
good_ids = (person_vs_pheno['obesitycat']!='Overweight').values & ((abundance>0).sum(1)>10).values

person_vs_pheno = person_vs_pheno[good_ids]

# Taxa aggregated by taxa category
person_variant_species = np.array(abundance.transpose().join(tax_table).groupby('Species').aggregate(sum).transpose())
person_variant_genus = np.array(abundance.transpose().join(tax_table).groupby('Genus').aggregate(sum).transpose())
person_variant_family = np.array(abundance.transpose().join(tax_table).groupby('Family').aggregate(sum).transpose())
person_variant_order = np.array(abundance.transpose().join(tax_table).groupby('Order').aggregate(sum).transpose())
person_variant_class = np.array(abundance.transpose().join(tax_table).groupby('Class').aggregate(sum).transpose())
person_variant_phylum = np.array(abundance.transpose().join(tax_table).groupby('Phylum').aggregate(sum).transpose())
person_variant_kingdom = np.array(abundance.transpose().join(tax_table).groupby('Kingdom').aggregate(sum).transpose())

person_variant_1 = np.array(np.load(PROJECT_DIR + '/results/person_biomarker/' + study  + 
                                   '_'  + str(msa_method) +  '_person_variant1_unique.npy'))
person_variant_2 = np.array(np.load(PROJECT_DIR + '/results/person_biomarker/' + study  + 
                                   '_'  + str(msa_method) +  '_person_variant2_unique.npy'))
person_variant_3 = np.array(np.load(PROJECT_DIR + '/results/person_biomarker/' + study  + 
                                   '_'  + str(msa_method) +  '_person_variant3_unique.npy'))

true_phenos = person_vs_pheno['pheno'].values



person_vs_pheno_dict = {i:j for i,j in zip(person_vs_pheno['family'], person_vs_pheno['pheno'])}
k = list(person_vs_pheno_dict.keys())
v = list(person_vs_pheno_dict.values())
f_FDR,ax_FDR = plt.subplots(figsize=(10,5))
f_FDR2,ax_FDR2 = plt.subplots(figsize=(10,5))
f_qq,ax_qq = plt.subplots(figsize=(10,5))

for person_variant, title in [(person_variant_species, 'species'), (person_variant_genus, 'genus'),  (person_variant_family, 'family'),  (person_variant_order, 'order'), (person_variant_class, 'class'),
    (person_variant_1, '1-SBB'),(person_variant_2, '2-SBB'),(person_variant_3, '3-SBB')]:
    print(title)
    
    # Compute Statistics
    person_variant = person_variant[good_ids]
    
    # Statistics
    U_data = np.zeros(np.shape(person_variant)[1])
    n1 = sum(true_phenos)
    N = len(true_phenos)
    print(sum(true_phenos))
    
    def computeMannU(i, phenos):
        if (i%100000==0): print(i)
        vals = person_variant[:, i]
        if sum(vals)==0:
            return FAKE_P_VAL
        _, pval = mannwhitneyu(vals[phenos], vals[~phenos])
        if np.isnan(pval):
            return FAKE_P_VAL
        else: return pval

    U_data = np.array(Parallel(n_jobs=num_cores-2)(delayed(computeMannU)(i, true_phenos)
                                                        for i in range(np.shape(person_variant)[1])))
    print("Done with U_Data")
    U_data_na = U_data[U_data!=FAKE_P_VAL]
    
    U_permute = np.zeros((np.shape(person_variant)[1], N_ITER))
    n_permute_gt_r = np.zeros(np.shape(person_variant)[1])
    for r in tqdm.tqdm(range(N_ITER)):
        np.random.shuffle(v)    
        new_pheno_dict = {i:j for i,j in zip(k,v)}
        new_phenos = np.array([new_pheno_dict[i] for i in person_vs_pheno['family']])
        while sum(new_phenos)!=sum(true_phenos): ## Make sure there are always the same # of obese vs lean
            np.random.shuffle(v)    
            new_pheno_dict = {i:j for i,j in zip(k,v)}
            new_phenos = np.array([new_pheno_dict[i] for i in person_vs_pheno['family']])
        U_permute[:,r] = Parallel(n_jobs=num_cores-2)(delayed(computeMannU)(i, new_phenos)
                                                        for i in range(np.shape(person_variant)[1]))
    U_permute = U_permute.flatten()
    U_permute_na = U_permute[U_permute!=FAKE_P_VAL]

    #U_threshs = np.linspace(np.percentile(U_data, 1),0,10000)
    #fdr = [np.mean(U_permute_na<=U_thresh)/np.mean(U_data_na<=U_thresh) for U_thresh in U_threshs]
    #num_hits = [np.sum(U_data_na<=U_thresh) for U_thresh in U_threshs]
    ##frac_hits = [np.mean(U_data_na<=U_thresh) for U_thresh in U_threshs]
    #fdr = np.array([max(fdr[f:]) for f in range(len(fdr))])
    #fdr, fdr_idx = np.unique(fdr, return_index=True)
    ##num_hits = np.array(num_hits)[fdr_idx]
    #frac_hits = np.array(frac_hits)[fdr_idx]
    
    #ax_FDR.plot(fdr, num_hits, 'x-', label=title)
    #ax_FDR2.plot(fdr, frac_hits, 'x-', label=title)

    #pvals = n_permute_gt_r/N_ITER
    #reject, adjusted_pvals, _, _ = multipletests(pvals, alpha=0.2, method='fdr_bh', is_sorted=False, returnsorted=False)
    
    #### Plot histogram of R values. ####
    #plt.figure()
    #plt.hist(U_data_na, alpha=.5, density=True, bins=np.linspace(0,1,100))
    ##plt.hist(U_permute_na, alpha=.5, density=True, bins=np.linspace(0,1,100))
    #plt.xlabel('P Value')
    #plt.ylabel('Freq')
    #plt.savefig(PLOT_DIR + 'yogurt_' + title + '_histogram.png')
    
    
    #quantiles_data = np.array([np.percentile(U_data, p) for p in np.linspace(0,100,1000)])
    ##quantiles_permute = np.array([np.percentile(U_permute, p) for p in np.linspace(0,100,1000)])
    #ax_qq.plot(quantiles_permute, quantiles_data, '-', label=title)
    #ax_qq.plot(quantiles_permute, quantiles_permute, 'k--')
    
    np.save(PROJECT_DIR + '/results/sig_values/obese_lean_twins_' + 'U_permute_' + title + '.npy',U_permute_na)
    np.save(PROJECT_DIR + '/results/sig_values/obese_lean_twins_' + 'U_data_' + title + '.npy',U_data_na)

    #i_sig = np.argmin(U_data)
    #print(min(pvals)*len(pvals))
    #df = pd.DataFrame([person_variant[:,i_sig], ['Obese' if t else 'Lean' for t in true_phenos]]).transpose()
    #df.columns = ['val', 'pheno']
    #df['val'] = df['val'].astype(float)
    
#ax_FDR.set_yscale('log')    
#ax_FDR.set_ylabel('# Significant SBBs')  # we already handled the x-label with ax1
#ax_FDR.set_xlabel('FDR')  # we already handled the x-label with ax1
#ax_FDR.set_xlim((0,.2))  # we already handled the x-label with ax1

#ax_FDR2.set_ylabel('Fraction Significant SBBs')  # we already handled the x-label with ax1
#ax_FDR2.set_xlabel('FDR')  # we already handled the x-label with ax1
#ax_FDR2.set_yscale('log')    
#ax_FDR2.set_xlim((0,.2))  # we already handled the x-label with ax1


##ax_qq.set_ylabel('True Quantile' )  # we already handled the x-label with ax1
#ax_qq.set_xlabel('Theoretical Quantile')  # we already handled the x-label with ax1
#f_qq.legend()
#f_qq.savefig(PLOT_DIR + 'obese_QQ.png')
#f_qq.show()

#f_FDR.legend()
#f_FDR.savefig(PLOT_DIR + 'obese_FDR_num.png')
#f_FDR.show()
#f_FDR2.legend()
#f_FDR2.savefig(PLOT_DIR + 'obese_FDR_frac.png')
#f_FDR2.show()
