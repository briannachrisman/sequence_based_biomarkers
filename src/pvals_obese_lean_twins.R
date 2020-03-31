source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/helpers.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/statistics.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/variant_combos.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/pipelines.R")
library(phyloseq)
library(msa)
library(RcppCNPy)
library(Matrix)

msa_alg = 'DECIPHER'
# Read in and format dataframes.
ps = readRDS("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/final/ps.rds")
MSA = readRDS(paste0("/oak/stanford/groups/dpwall//users/briannac/phyloWAS/data/obese_lean_twins/msa_", msa_alg, ".rds"))
person_biomarker_dir = paste0("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/obese_lean_twins/person_biomarker/", msa_alg, "/")
save_dir = paste0("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/obese_lean_twins/p_vals/", msa_alg, "_")


# Transform sample counts and merge appropriately.
print("Transforming Sample Counts...")
ps = transform_sample_counts(ps, function(x){x/sum(x)})
options(warn=-1)
ps_sub = merge_samples(ps, 'host_subject_id')  # Change to whatever column contains host info!
options(warn=1)
sample_names(ps_sub) = sample_names(ps)[!duplicated(sample_data(ps)$host_subject_id)] # Change to whatever column contains host info!
sample_data(ps_sub) = sample_data(ps)[!duplicated(sample_data(ps)$host_subject_id),] # Change to whatever column contains host info!
ps = ps_sub

print("Transforming Sample Counts...")
ps = transform_sample_counts(ps, function(x){x/sum(x)})
seq_df = t(as.data.frame(sapply(as.character(MSA), function(x){strsplit(x, "")})))
rownames(seq_df) = taxa_names(ps_sub)
seq_df = seq_df[taxa_names(ps),]
                                            
### Perform phylo-was.
n_permutations = 100
phenotypes = sample_data(ps)$obesitycat == 'Obese' ## TO CHANGE: Change to whatever column has relevant phenotype!

bases = c('A','T','C','G')
n_cores = 15
idx = sample_data(ps)$family ##### TO CHANGE: Change to whatever column has family info! (Or if no family info, then just do host.)
paired = F

##### Combos of 1. ######
n_df = 1

print("PhyloWAS on combos of 1...")
person_biomarker = npyLoad(paste0(person_biomarker_dir, 'person_variant', n_df, '_condensed.npy'))
var_names = read.table(paste0(person_biomarker_dir, 'biomarkers', n_df, '.txt'))
              
              
exports = RunOnBiomarkers(person_variant = person_biomarker, var_names = var_names, phenotypes = phenotypes, group_idx = idx,
                              n_df = n_df, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
exports$seq_df = seq_df
saveRDS(exports, paste0(save_dir, n_df, '.rds'))


##### Combos of 2. ######
n_df = 2
print(paste0("PhyloWAS on combos of ", n_df, "..."))
person_biomarker = npyLoad(paste0(person_biomarker_dir, 'person_variant', n_df, '_condensed.npy'))
var_names = read.table(paste0(person_biomarker_dir, 'biomarkers', n_df, '.txt'))
exports = RunOnBiomarkers(person_variant = person_biomarker, var_names = var_names, phenotypes = phenotypes, group_idx = idx,
                              n_df = n_df, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
exports$seq_df = seq_df
saveRDS(exports, paste0(save_dir, n_df, '.rds'))


##### Combos of 3. ######
#n_df = 3
#print(paste0("PhyloWAS on combos of ", n_df, "..."))
#person_biomarker = npyLoad(paste0(person_biomarker_dir, n_df, '.npy'))
#var_names = read.table(paste0(person_biomarker_dir, n_df, '_names.txt'))
#exports = RunOnBiomarkers(person_variant = person_biomarker, var_names = var_names, phenotypes = phenotypes, group_idx = idx,
#                              n_df = n_df, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
#exports$ps = ps
#exports$seq_df = seq_df
#saveRDS(exports, paste0(save_dir, n_df, '.rds'))