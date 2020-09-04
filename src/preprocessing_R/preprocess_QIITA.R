source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/preprocessing.R")

# Cardiometabolic study from Columbia.
sample_ids = read.csv("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_colombia/split/sample_ids.txt", header=F)
sample_ids = unique(as.vector(as.character(sample_ids[["V1"]])))
out_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_colombia/split/"
save_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_colombia/preprocessing_out/"
final_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_colombia/final/"
fastq_all_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_colombia/fastq/Demultiplexed/62374/seqs.fastq"
mapping_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/cardiometabolic_colombia/split/mapping.csv"

# Chron's.
#sample_ids = read.csv("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/split/sample_ids.txt", header=F)
#sample_ids = unique(as.vector(as.character(sample_ids[["V1"]])))
#out_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/split/"
#save_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/preprocessing_out/"
#final_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/final/"
#fastq_all_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/fastq/Demultiplexed/29803/seqs.fastq"
#mapping_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/crohns/split/mapping.csv"



# Obese Lean Twins
#sample_ids = read.csv("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/split/sample_ids.txt", header=F)
#sample_ids = unique(as.vector(as.character(sample_ids[["V1"]])))
#out_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/split/"
#save_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/preprocessing_out/"
#final_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/final/"
#fastq_all_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/fastq/Demultiplexed/6821/seqs.fastq"
#mapping_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/obese_lean_twins/split/mapping.csv"

# Type 1 Diabetes Study.
#sample_ids = read.csv("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes/split/sample_ids.txt", header=F)
#sample_ids = sample_ids[sample_ids!='11129.NT1D033.JJK.St']
#out_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes/split/"
#save_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes/preprocessing_out/"
#final_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes/final/"
#fastq_all_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes/fastq/Demultiplexed/31386/seqs.fastq"
#mapping_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes/split/mapping.csv"

# Type 1 Diabetes Study -- chopped 150
#sample_ids = read.csv("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes_150/split/sample_ids.txt", header=F)
#sample_ids = unique(as.vector(as.character(sample_ids[["V1"]])))
#sample_ids = sample_ids[sample_ids!='11129.NT1D033.JJK.St']
#out_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes_150/split/"
#save_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes_150/preprocessing_out/"
#final_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes_150/final/"
#fastq_all_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes_150/fastq/Demultiplexed/31387/seqs.fastq"
#mapping_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/diabetes_150/split/mapping.csv"

# IBD Twins Study
#sample_ids = read.csv("/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/jansson_twins_ibd/split/sample_ids.txt", header=F)
#sample_ids = unique(as.vector(as.character(sample_ids[["V1"]])))
#out_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/jansson_twins_ibd/split/"
#save_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/jansson_twins_ibd/preprocessing_out/"
#final_dir = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/jansson_twins_ibd/final/"
#fastq_all_file = '/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/jansson_twins_ibd/fastq/Demultiplexed/6485/seqs.fastq'
#mapping_file = "/oak/stanford/groups/dpwall/microbiome/ncbi-datasets/jansson_twins_ibd/split/mapping.csv"

trimLeft = 0
truncLen = 150

dadaPipeline(sample_ids, out_dir, save_dir, final_dir, fastq_all_file, mapping_file, trimLeft, truncLen, minBoot=50,
         ref_species_file = "/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/ref_genome/rdp_species_assignment_14.fa.gz",
             ref_taxonomy_file = "/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/ref_genome/rdp_train_set_16.fa.gz")