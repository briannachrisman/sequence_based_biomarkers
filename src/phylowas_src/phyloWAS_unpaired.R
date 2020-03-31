source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/helpers.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/statistics.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/variant_combos.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/pipelines.R")
library(phyloseq)
library(msa)

# Read in and format dataframes.
ps = readRDS("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/final/ps.rds")
MSA = readRDS("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/data/crohns/msa_DECIPHER.rds")
save_dir = "/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/crohns/phyloWAS_DECIPHER"

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
#ps = subset_taxa(ps, taxa_names(ps) %in% taxa_names(ps)[order(colSums(otu_table(ps)), decreasing=T)[1:20]]) ## Subset taxa just for testing.
ps = transform_sample_counts(ps, function(x){x/sum(x)})
seq_df = t(as.data.frame(sapply(as.character(MSA), function(x){strsplit(x, "")})))
rownames(seq_df) = taxa_names(ps_sub)
seq_df = seq_df[taxa_names(ps),]
consensus_mat = consensusMatrix(MSA)
consensus_seq = apply(consensus_mat, 2, function(x) rownames(consensus_mat)[which.max(x)])
isConsensus = t(apply(seq_df, 1, function(x) (x==consensus_seq)))
seq_df[isConsensus] = 'X'
                      
                      
### Perform phylo-was.
n_permutations = 100
phenotypes = sample_data(ps)$diseasesubtype == 'UC' ## TO CHANGE: Change to whatever column has relevant phenotype!
bases = c('A','T','C','G')
n_cores = 15
idx = sample_data(ps)$host_subject_id ##### TO CHANGE: Change to whatever column has family info! (Or if no family info, then just do host.)
paired = F

# Single taxa.
print("PhyloWAS on single taxa...")
person_taxa = as.sparseMatrix(otu_table(ps), use_pointers = T)
taxa_variant_combos = sparseMatrix(1:ncol(person_taxa), 1:ncol(person_taxa), dims=c(ncol(person_taxa), ncol(person_taxa)))
exports = RunSingleVariantMWAS(person_taxa = person_taxa, taxa_variant_combos = taxa_variant_combos,
                               phenotypes = phenotypes, idx = idx, paired = paired,
                               n_permutations=n_permutations, n_cores = n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_ASVs.rds'))

# Tax level.
for (taxa_level in c("Phylum", "Order", "Class", "Family", "Genus", "Species")) {
    print(paste0("PhyloWAS on", taxa_level, "..."))
    ps_tax_level = tax_glom(ps, taxa_level)
    person_taxa = as.sparseMatrix(otu_table(ps_tax_level), use_pointers = T)
    taxa_variant_combos = sparseMatrix(1:ncol(person_taxa), 1:ncol(person_taxa), dims=c(ncol(person_taxa), ncol(person_taxa)))
    exports = RunSingleVariantMWAS(person_taxa = person_taxa, taxa_variant_combos = taxa_variant_combos,
                               phenotypes = phenotypes, idx = idx, paired = paired,
                               n_permutations=n_permutations, n_cores = n_cores)
    exports$ps = ps_tax_level
    saveRDS(exports, paste0(save_dir, '_', taxa_level, '.rds'))
}
person_taxa = as.sparseMatrix(otu_table(ps), use_pointers = T)

###### Combos of 2. #####
print("PhyloWAS on combos of 2 basepairs...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, idx = idx,
                              n_df = 2, paired = paired, bases = bases, n_permutations=n_permutations, n_cores = n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_2.rds'))


##### Combos of 1. ######
print("PhyloWAS on combos of 1 basepairs...")
exports = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, idx = idx,
                              n_df = 1, paired = paired, bases = bases, n_permutations=n_permutations, n_cores=n_cores)
exports$ps = ps
saveRDS(exports, paste0(save_dir, '_combos_of_1.rds'))
print("Done --- success!!! :)")
