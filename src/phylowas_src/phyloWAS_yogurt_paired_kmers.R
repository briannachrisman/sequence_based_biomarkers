source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/helpers.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/statistics.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/variant_combos.R")
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/pipelines.R")
library(phyloseq)
library(DECIPHER)

n_permutations = 100
n_cores = 15
wordSize = 5
save_dir = paste0("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/results/yogurt/phyloWAS_paired_", wordSize, 'mrs.rds')

# Read in and format dataframes.
ps_seq = readRDS("/oak/stanford/groups/dpwall//users/briannac/Microbiome/data/project_yogurt/phyloseq_20190611.rds") # These are numbered weirdly. The code below fixes it (debugged 10/6/2019).
ps_og = readRDS("/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/project_yogurt/phyloseq_20190909.rds")

samples = sample_data(ps_seq)
taxa_names(ps_seq) = paste0('ASV_', 1:ntaxa(ps_seq))


ps = ps_og
ps = transform_sample_counts(ps, function(x){x/sum(x)})
options(warn=-1)
ps_sub = merge_samples(ps, 'host')
options(warn=1)
sample_names(ps_sub) = sample_names(ps)[!duplicated(sample_data(ps)$host)]
sample_data(ps_sub) = sample_data(ps)[!duplicated(sample_data(ps)$host),]
ps = ps_sub
#ps = subset_samples(ps, sample_data(ps)$age>24) # Only keep children over 24 months.
x=table(paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint))
ps = subset_samples(ps, paste0(sample_data(ps)$family, '_', sample_data(ps)$timepoint) %in% names(x)[which(x==2)]) # Only keep PAIRED samples.
ps = subset_taxa(ps, colSums(otu_table(ps))>0) # Only keep PAIRED samples.
ps = subset_taxa(ps, taxa_names(ps) %in% c('ASV_403', 'ASV_1', taxa_names(ps)[order(colSums(otu_table(ps)), decreasing=T)[1:20]])) ## Subset taxa just for testing.  # ASV_1 and ASV_403 have super similar seqs.
ps = transform_sample_counts(ps, function(x){x/sum(x)})

# Perform Multiple Sequence Alignment.
ps_seq = subset_taxa(ps_seq, taxa_names(ps_seq) %in% taxa_names(ps))
seqs = tax_table(ps_seq)
seqs = seqs[taxa_names(ps_seq),]
seq_set = DNAStringSet(seqs[,'Sequence'])
kmers = .Call("enumerateSequence", seq_set, wordSize, PACKAGE="DECIPHER")
kmer_melt = melt(kmers)
kmer_mat = sparseMatrix(i = kmer_melt$L1, j = kmer_melt$value+1, dims = c(length(kmers), max(kmer_melt$value)+1))
print(paste0(sum(apply(kmer_mat, 2, sum)>1), " kmers."))
kmer_mat = kmer_mat[, apply(kmer_mat, 2, sum)>1]
kmer_mat = kmer_mat[,!duplicated.sparseMatrix(kmer_mat, 2)]
taxa_variant_combos = kmer_mat                      
                      
### Perform phylo-was.
phenotypes = sample_data(ps)$phenotype == 'A'
idx = sample_data(ps)$family
paired = T

# K-mers.
person_taxa = as.sparseMatrix(otu_table(ps), use_pointers = T)
taxa_variant_combos = sparseMatrix(1:ncol(person_taxa), 1:ncol(person_taxa), dims=c(ncol(person_taxa), ncol(person_taxa)))
exports = RunSingleVariantMWAS(person_taxa = person_taxa, taxa_variant_combos = taxa_variant_combos,
                               phenotypes = phenotypes, idx = idx, paired = paired,
                               n_permutations=n_permutations, n_cores = n_cores)
exports$ps = ps
saveRDS(exports, save_dir)