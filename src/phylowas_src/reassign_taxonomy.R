library(BiocParallel)
library(doParallel)
register(MulticoreParam(workers=15))
registerDoParallel(cores=15)
library(phyloseq)
library(dada2)
seq_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/final/seqs.rds"
ps_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/final/ps.rds"
new_ps_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/final/ps_with_taxa.rds"
ref_taxonomy_file = "/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/ref_genome/GTDB_bac-arc_ssu_r86.fa.gz"
seqs = readRDS(seq_file)
seq_names = seqs[,'name']
colnames(seqs) = c('sequence', 'abundance')
ps = readRDS(ps_file)
print("Assigning Taxonomy...")
taxa = assignTaxonomy(seqs, ref_taxonomy_file, multithread = T, verbose=T)
row.names(taxa) = seq_names
tax_table(ps) = taxa
saveRDS(ps, new_ps_file)
