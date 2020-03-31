library(phyloseq)
library(Biostrings)
library(msa)
library(DECIPHER)


# For Autism dataset.
#save_dir = '/oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/yogurt/'
#ps = readRDS('/oak/stanford/groups/dpwall//users/briannac/Microbiome/data/project_yogurt/phyloseq_20190829.rds')
#ps_og = readRDS('/oak/stanford/groups/dpwall//users/briannac/Microbiome/data/project_yogurt/phyloseq_20190611.rds') # These are numbered weirdly. The code below fixes it #(debugged 10/6/2019).
#samples = sample_data(ps_og)
#taxa_names(ps_og) = paste0('ASV_', 1:ntaxa(ps_og))
#rownames(samples) = sample_names(ps)
#sample_data(ps) = cbind(sample_data(ps), samples)
#seqs = tax_table(ps_og)

# For obese/lean twins dataset.
save_dir = '/oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/crohns/'
ps = readRDS('/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/final/ps.rds')
seqs = readRDS('/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/final/seqs.rds')

# For cardiometabolic.
#save_dir = '/oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/cardiometabolic_columbia'
#ps = readRDS('/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_columbia/final/ps.rds')
#seqs = readRDS('/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_columbia/final/seqs.rds')


# Perform Multiple Sequence Alignment.
seq_set = DNAStringSet(seqs[,'Sequence'])
print("Performing DECIPHER MSA....")
aligned = AlignSeqs(seq_set, anchor=NA,verbose=T)
saveRDS(aligned, paste0(save_dir, 'msa_DECIPHER.rds'))

print("Performing ClustalW MSA....")
aligned = msa(seq_set, method="ClustalW")
saveRDS(aligned, paste0(save_dir, 'msa_ClustalW.rds'))

print("Performing ClustalOmega MSA....")
aligned = msa(seq_set, method="ClustalOmega")
saveRDS(aligned, paste0(save_dir, 'msa_ClustalOmega.rds'))

print("Performing Muscle MSA....")
aligned = msa(seq_set, method="Muscle")
saveRDS(aligned, paste0(save_dir, 'msa_Muscle.rds'))
print("Done --- success! :)")