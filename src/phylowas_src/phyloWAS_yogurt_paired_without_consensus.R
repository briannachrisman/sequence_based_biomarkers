source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/helpers.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/statistics.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/variant_combos.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/pipelines.R")
library(phyloseq)
library(msa)

# Read in and format dataframes.
ps_og = readRDS("/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/project_yogurt/phyloseq_20190909.rds")
MSA = readRDS("/oak/stanford/groups/dpwall//users/briannac/phyloWAS/data/yogurt/msa_DECIPHER.rds")
save_dir = "/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/yogurt/phyloWAS_DECIPHER_without_consensus"

print("Reading in data...")
ps = ps_og
ps = transform_sample_counts(ps, function(x){x/sum(x)})
options(warn=-1)
ps_sub = merge_samples(ps, 'host')
options(warn=1)
sample_names(ps_sub) = sample_names(ps)[!duplicated(sample_data(ps)$host)]
sample_data(ps_sub) = sample_data(ps)[!duplicated(sample_data(ps)$host),]
ps = ps_sub

print("Cleaning up...")
#ps = subset_samples(ps, sample_data(ps)$age>24) # Only keep children over 24 months.
x=table(paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint))
ps = subset_samples(ps, paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint) %in% names(x)[which(x==2)]) # Only keep PAIRED samples.
ps = subset_taxa(ps, colSums(otu_table(ps))>0) # Only keep PAIRED samples.

#ps = subset_taxa(ps, taxa_names(ps) %in% c('ASV_403', 'ASV_1', taxa_names(ps)[order(colSums(otu_table(ps)), decreasing=T)[1:20]])) ## Subset taxa just for testing.  # ASV_1 and ASV_403 have super similar seqs.

message("Dealing with msa...")
ps = transform_sample_counts(ps, function(x){x/sum(x)})
seq_df = t(as.data.frame(sapply(as.character(MSA), function(x){strsplit(x, "")})))
rownames(seq_df) = taxa_names(ps_og)
seq_df = seq_df[taxa_names(ps),]
consensus_mat = consensusMatrix(MSA)
consensus_seq = apply(consensus_mat, 2, function(x) rownames(consensus_mat)[which.max(x)])
isConsensus = t(apply(seq_df, 1, function(x) (x==consensus_seq)))
seq_df[isConsensus] = 'X'
                      
                      
### Perform phylo-was.
message("PhyloWASing...")

n_permutations = 100
phenotypes = sample_data(ps)$phenotype == 'A'
bases = c('A','T','C','G')
n_cores = 19
idx = sample_data(ps)$family
paired = T

# Single taxa.
message("Converting to sparse representtaions...")
person_taxa = as(otu_table(ps), "sparseMatrix")    
#taxa_variant_combos = sparseMatrix(1:ncol(person_taxa), 1:ncol(person_taxa), dims=c(ncol(person_taxa), ncol(person_taxa)))
#message("PhyloWAS on ASVs...")
#exports = RunSingleVariantMWAS(person_taxa = person_taxa, taxa_variant_combos = taxa_variant_combos,
#                               phenotypes = phenotypes, group_idx = idx, paired = paired,
#                               n_permutations=n_permutations, n_cores = n_cores)
#exports$ps = ps
#saveRDS(exports, paste0(save_dir, '_ASVs.rds'))

# Tax level.

#for (taxa_level in c("Phylum", "Order", "Class", "Family", "Genus", "Species")) {
#    print(paste0("PhyloWAS on ", taxa_level, "..."))
#    ps_tax_level = tax_glom(ps, taxa_level)
#    person_taxa = as.sparseMatrix(otu_table(ps_tax_level), use_pointers = T)
#    taxa_variant_combos = sparseMatrix(1:ncol(person_taxa), 1:ncol(person_taxa), dims=c(ncol(person_taxa), ncol(person_taxa)))
#    exports = RunSingleVariantMWAS(person_taxa = person_taxa, taxa_variant_combos = taxa_variant_combos,
#                               phenotypes = phenotypes, group_idx = idx, paired = paired,
#                               n_permutations=n_permutations, n_cores = n_cores)
#    exports$ps = ps_tax_level
#    saveRDS(exports, paste0(save_dir, '_', taxa_level, '.rds'))
#}
#person_taxa = as.sparseMatrix(otu_table(ps), use_pointers = T)

##### Combos of 1. ######
print("PhyloWAS on combos of 1...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 1, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_1.rds'))
                      
                      
###### Combos of 2. #####
print("PhyloWAS on combos of 2...")
#
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 2, paired = paired, bases = bases, n_permutations=n_permutations, n_cores = n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_2.rds'))

                      
### Combos of 3 ######
message("PhyloWAS on combos of 3...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 3, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_3.rds'))
                      
                      
print("PhyloWAS on combos of 4...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 4, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_4.rds'))
                      
                      
print("PhyloWAS on combos of 5...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 5, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_5.rds'))
                      
                      
print("PhyloWAS on combos of 6...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 6, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_6.rds'))
                      
                      
                      
print("PhyloWAS on combos of 7...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 7, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
#saveRDS(exports, paste0(save_dir, '_combos_of_7.rds'))
                      
                      
print("PhyloWAS on combos of 8...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, group_idx = idx,
                              n_df = 5, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_8.rds'))