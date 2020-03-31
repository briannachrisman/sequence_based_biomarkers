source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/helpers.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/statistics.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/variant_combos.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/pipelines.R")
library(phyloseq)
library(msa)
library(RcppCNPy)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
msa_alg = args[1]
# Read in and format dataframes.
MSA = readRDS(paste0("/oak/stanford/groups/dpwall//users/briannac/phyloWAS/data/yogurt/msa_", msa_alg, ".rds"))
person_biomarker_dir = paste0("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/yogurt/person_biomarker/", msa_alg, "/")
save_dir = paste0("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/yogurt/p_vals/", msa_alg, "_discrete_times")
ps_og = readRDS("/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/project_yogurt/phyloseq_20190909.rds")

print("Reading in data...")
ps = ps_og
ps = transform_sample_counts(ps, function(x){x/sum(x)})
options(warn=-1)
#ps_sub = merge_samples(ps, 'host')
options(warn=1)
#sample_names(ps_sub) = sample_names(ps)[!duplicated(sample_data(ps)$host)]
#sample_data(ps_sub) = sample_data(ps)[!duplicated(sample_data(ps)$host),]
#ps = ps_sub

print("Cleaning up...")
x=table(paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint))
ps = subset_samples(ps, paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint) %in% names(x)[which(x==2)]) # Only keep PAIRED samples.
ps = subset_taxa(ps, colSums(otu_table(ps))>0) # Only keep PAIRED samples.

print("Dealing with MSA...")
seq_df = t(as.data.frame(sapply(as.character(MSA), function(x){strsplit(x, "")})))
rownames(seq_df) = taxa_names(ps_og)
seq_df = seq_df[taxa_names(ps),]

### Perform phylo-was.
message("PhyloWASing...")

n_permutations = 100
phenotypes = sample_data(ps)$phenotype == 'A'
bases = c('A','T','C','G')
n_cores = 19
idx = paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint)
paired = T

##### Combos of 1. ######
#n_df = 1

#print("PhyloWAS on combos of 1...")
#person_biomarker = npyLoad(paste0(person_biomarker_dir, 'person_variant', n_df, '_condensed.npy'))
#var_names = read.table(paste0(person_biomarker_dir, 'biomarkers', n_df, '.txt'))
#exports = RunOnBiomarkers(person_variant = person_biomarker, var_names = var_names, phenotypes = phenotypes, group_idx = idx,
#                              n_df = n_df, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
#exports$ps = ps
#exports$seq_df = seq_df
#saveRDS(exports, paste0(save_dir, n_df, '.rds'))


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
print('success!!')

##### Combos of 3. ######
#n_df = 3
#print(paste0("PhyloWAS on combos of ", n_df, "..."))
#person_biomarker = npyLoad(paste0(person_biomarker_dir, 'person_variant', n_df, '_condensed.npy'))
#var_names = read.table(paste0(person_biomarker_dir, 'biomarkers', n_df, '.txt'))
#exports = RunOnBiomarkers(person_variant = person_biomarker, var_names = var_names, phenotypes = phenotypes, group_idx = idx,
#                              n_df = n_df, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
#exports$ps = ps
#exports$seq_df = seq_df
#saveRDS(exports, paste0(save_dir, n_df, '.rds'))
#print('success!!')