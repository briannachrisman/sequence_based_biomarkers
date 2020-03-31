library(dada2)
library(phyloseq)
library(BiocParallel)
library(doParallel)
register(MulticoreParam(workers=19))
registerDoParallel(cores=19)
library(ShortRead)
library(data.table)
library(plyr)

splitFastq = function(fastq_all, ids, sample_id ,out_dir) {
    print(paste0(out_dir, sample_id, '.fastq.gz'))
    idx = grep(paste0(sample_id, '_'), ids)
    reads = fastq_all[idx]
    writeFastq(object = reads , file = paste0(out_dir, sample_id, '.fastq.gz'), compress = T)
}

dadaPipeline = function(sample_ids, out_dir, save_dir, final_dir, fastq_all_file, mapping_file, trimLeft, truncLen, minBoot=50, max_seqs_in_group = 6,
                        ref_species_file = "/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/ref_genome/rdp_species_assignment_14.fa.gz",
                        ref_taxonomy_file = "/oak/stanford/groups/dpwall/users/briannac/Microbiome/data/ref_genome/rdp_train_set_16.fa.gz") {
    # Split FastQ.
    print("Reading fastq file...")
    fastq_all = readFastq(fastq_all_file)
    ids = ShortRead::id(fastq_all)
    groups = c(seq(0,length(sample_ids), max_seqs_in_group+1), length(sample_ids))
    sample_ids_groups = sapply(2:length(groups), function(i) {sample_ids[(groups[i-1]+1):groups[i]]})
    i_group = 0
    full_seq_table = data.frame()
    sample_ids_list = c()
    for (sample_ids_group in sample_ids_groups) { ### Change for full pipeline (set to [1:2] for testing)
        sample_ids_list = c(sample_ids_group, sample_ids_list)
        print(sample_ids_group)
        print("Splitting fastq...")
        print(paste0("On group ", i_group))
        n = length(sample_ids_group)
        print(n)
        #print("Splitting fastq file...")
        #foreach(i=1:n) %dopar% splitFastq(fastq_all, ids, sample_ids_group[i], out_dir)
        i_group = i_group + 1
        fastq_files = as.vector(sapply(sample_ids_group, function(x) {paste0(out_dir, x, '.fastq.gz')}))
        filt_files =as.vector(sapply(sample_ids_group, function(x) {paste0(out_dir, x, '_filt.fastq.gz')}))
        
        # Quality Profile.
        print("Plotting quality profile...")
        #plotQualityProfile(fastq_files)
        #png(filename = paste0(save_dir,'_', i_group, '_quality_profile.png'))
    
        # Filter/Trim
        print("Filtering/trimming...")
        #filt_out = filterAndTrim(fastq_files, filt_files, trimLeft = trimLeft, truncLen = truncLen, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

        # Errors 
        print("Learning errors...")
        err = learnErrors(filt_files, multithread=T, verbose=T)
        plotErrors(err, nominalQ=TRUE)
        png(filename = paste0(save_dir, '_', i_group, '_errors.png'))

        # Derep / dada
        print("Derepping...")
        derep = derepFastq(filt_files, verbose=TRUE)
        print("Running dada...")
        dada_reads = dada(derep, err=err, multithread=TRUE, verbose = T)
        seq_table = makeSequenceTable(dada_reads)
        
        row_names = c(row.names(seq_table), row.names(full_seq_table))
        full_seq_table = join(data.frame(seq_table), full_seq_table, type='full')
        row.names(full_seq_table) = row_names
    }
    full_seq_table[is.na(full_seq_table)]=0
    saveRDS(full_seq_table, paste0(final_dir, 'seq_table.rds'))
    #full_seq_table = readRDS(full_seq_table)
    #colnames(full_seq_table) = c('sequence', 'abundance')

    print("Remove Bimeras...")
    no_chimeras = removeBimeraDenovo(as.matrix(full_seq_table), method="consensus", multithread=TRUE, verbose=TRUE)
    row.names(no_chimeras) = sample_ids_list
    

    # Read in mapping file.
    print("Dealing with mapping file...")
    mapping = read.csv(mapping_file, sep='\t')
    row.names(mapping) = mapping[,1]
    mapping = mapping[row.names(no_chimeras),]

    # Deal with sequence file.
    print("Dealing with sequence file...")
    seqs = data.frame(Sequence=colnames(no_chimeras))
    seqs$name = sapply(1:nrow(seqs), function(x) {paste0('ASV_', x)})
    row.names(seqs) = seqs$Sequence
    
    # Assign taxonomy.
    print("Assigning Taxonomy...")
    taxa = assignTaxonomy(no_chimeras, ref_taxonomy_file, multithread = T, minBoot=minBoot)
    taxa = addSpecies(taxa, ref_species_file)
    
    # Turn into phyloseq object + save.
    print("Converting to phyloseq and saving...")

    ps = phyloseq(otu_table(no_chimeras, taxa_are_rows=FALSE), 
               sample_data(mapping), tax_table(taxa))
    taxa_names(ps) = seqs$name
    saveRDS(ps, paste0(final_dir, 'ps.rds'))
    saveRDS(seqs, paste0(final_dir, 'seqs.rds'))
    print("Done -- Success!!! :)")
}