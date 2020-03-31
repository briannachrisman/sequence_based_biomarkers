library(reshape)
library(doParallel)
library(BiocParallel)
RunVariantComboMWAS = function(person_taxa, seq_df, phenotypes, group_idx, n_df, paired=T, bases=c(), n_permutations=10, n_cores=15) {
    register(MulticoreParam(workers=n_cores))
    registerDoParallel(cores=n_cores, )

    # Compute Person Variant Combos.
    message("Computing taxa_variants")
    t = Sys.time()
    t = SetTime(t)
    message("Computing taxa_variants -- converting to sparse matrix.")
    taxa_variant = TaxaVariantSparse(seq_df, variant_vals = bases)
    t = SetTime(t)

    message("Computing taxa_variants -- trimming.")
    to_keep = apply(taxa_variant>0, 2, sum)>0
    taxa_variant = taxa_variant[,to_keep]
    to_keep = !duplicated.sparseMatrix(taxa_variant, 2, has_vals=F)
    taxa_variant = taxa_variant[,to_keep]
    t = SetTime(t)

    #### Get sparse matrix represenation of taxa vs variants matrix.  ###
    message("Computing taxa_variants")
    # Table so we can easily switch between N-d ids to 1d id.
    #idx_lookup_table = array(1:(ncol(taxa_variant)**n_df), rep(ncol(taxa_variant), n_df))
    #message(paste0("\nSize of idx_lookup_table: ", object.size(idx_lookup_table)/1e9, " GB"))
    person_variants = sparseMatrix(i = c(1), j = c(1) , dims = c(ncol(taxa_variant)**n_df,nrow(person_taxa)))
    person_variants[1,1] = 0
    for (i_taxa in 1:nrow(taxa_variant)) {
        person_variants = person_variants + (t(TaxaVariantCombos(
        taxa_variant_sparse = taxa_variant, i_taxa = i_taxa, n_df = n_df)) %*% as.sparseVector(person_taxa[,i_taxa], T))
    }
#    person_variants =  (i_taxa = 1:nrow(taxa_variant), .combine='+') %dopar% (t(TaxaVariantCombos(
#        taxa_variant_sparse = taxa_variant, i_taxa = i_taxa, n_df = n_df)) %*% #as.sparseVector(person_taxa[,i_taxa], T)) 
    person_variants = t(person_variants)    
    var_ids=1:ncol(taxa_variant)
    t = SetTime(t)
    

    message("Getting variant names..")
    ### Getting names of variants. ###
    for (i in 1:n_df) {
        if (i==1) {df = data.frame(var_ids)}
        else {df = cbind(df, var_ids)}
    }
    idx = expand.grid(df)
    if (n_df >= 2) {

        for (i in 2:n_df) {
            person_variants = person_variants[,idx[,i-1] < idx[,i]]
            idx = idx[idx[,i-1] < idx[,i],]
        }
    }
    
    ## Trimming ##
    message(paste0("\nSize of person_variants: ", object.size(person_variants)/1e9, " GB"))
    message("\ntrimming taxa variants")
    message(paste0("\nTotal num of possible variants: ", ncol(person_variants)))

    # Get rid of biomarkers that don't exist in taxa set.
    to_keep = unique(summary(person_variants)$j)
    message(paste0("\nNon-zero variants: ", length(to_keep)))
    person_variants = person_variants[,to_keep]
    idx = as.matrix(idx[to_keep,])
    
    # Get rid of biomarkers that are perfectly in LD.
    to_keep = which(!duplicated.sparseMatrix(person_variants, 2, has_vals=T))
    message(paste0("\nNon-duplicate variants: ", length(to_keep)))
    person_variants = person_variants[,to_keep]
    idx = as.matrix(idx[to_keep,])
    
    message("\nConverting var_names_to_string")
    var_names = apply(idx, 1, function(x) {(paste(colnames(taxa_variant)[x], sep = '', collapse = "_"))})
    t = SetTime(t)

    # Compute p-values.
    t = SetTime(t)
    message("\nComputing p-values")
    if (paired) {
        permuted_pairs = PermuteDataset(ids = group_idx[which(phenotypes)], n_permute_iter = n_permutations)
        p_all = data.frame(
            foreach(i_group=1:ncol(person_variants))
            %dopar% ComputePairedP(person_variants, phenotypes = phenotypes, shuff_switch_sign = permuted_pairs, i_group = i_group))
        }
    else {
        permuted_phenos = PermuteDataset(ids = group_idx, n_permute_iter = n_permutations, prob=c(mean(phenotypes==T), mean(phenotypes==F)))
        p_all = data.frame(
            foreach(i_group=1:ncol(person_variants))
            %dopar% ComputeUnpairedP(person_variants, phenotypes = phenotypes, shuff_phenos = permuted_phenos, i_group = i_group))
    }
    
    # Set up exports.
    p = melt(p_all[1,])$value # P-values for the true data.
    p_shuff = p_all[-1,] # P-values for our shuffled data.
    t = SetTime(t=NA)
    message("Exporting Data.")
    export = list()
    export$p = p
    export$p_shuff = p_shuff
    export$taxa_variant = taxa_variant
    #export$taxa_variant_combos = taxa_variant_combos
    export$var_names = var_names
    export$person_variants = person_variants
    export$seq_df = seq_df

    if (paired) { export$permuted_pairs = permuted_pairs }
    else {export$permuted_phenos = permuted_phenos }
    export$person_taxa = person_taxa
    message("Success.")
    return(export)
}


RunSingleVariantMWAS = function(person_taxa, taxa_variants, phenotypes, group_idx, n_df, paired=T, bases=c(), n_permutations=10, n_cores=15) {
    taxa_variants = 
    register(MulticoreParam(workers=n_cores))
    registerDoParallel(cores=n_cores, )
    message("Computing person_variant.")
    person_variants = person_taxa %*% taxa_variants
    # Compute p-values.
    message("Computing p-values")
    if (paired) {
        permuted_pairs = PermuteDataset(ids = group_idx[which(phenotypes)], n_permute_iter = n_permutations)
        p_all = data.frame(
            foreach(i_group=1:ncol(person_variants))
            %dopar% ComputePairedP(person_variants, phenotypes = phenotypes, shuff_switch_sign = permuted_pairs, i_group = i_group))
        }
    else {
        permuted_phenos = PermuteDataset(ids = group_idx, n_permute_iter = n_permutations, prob=c(mean(phenotypes==T), mean(phenotypes==F)))
        p_all = data.frame(
            foreach(i_group=1:ncol(person_variants))
            %dopar% ComputeUnpairedP(person_variants, phenotypes = phenotypes, shuff_phenos = permuted_phenos, i_group = i_group))
    }
    
    # Set up exports.
    p = melt(p_all[1,])$value # P-values for the true data.
    p_shuff = p_all[-1,] # P-values for our shuffled data.
    t = SetTime(t=NA)
    message("Exporting Data.")
    export = list()
    export$p = p
    export$p_shuff = p_shuff
    export$taxa_variants = taxa_variants
    export$person_variants = person_variants
    if (paired) { export$permuted_pairs = permuted_pairs }
    else {export$permuted_phenos = permuted_phenos }
    export$person_taxa = person_taxa
    return(export)
    message("Success.")
}


RunOnBiomarkers = function(person_variants, var_names, phenotypes, group_idx, n_df, paired=T, bases=c(), n_permutations=10, n_cores=15) {
    to_keep = !duplicated(person_variants, MARGIN=2)
    person_variants = person_variants[,to_keep]
    var_names = var_names[to_keep,]
    
    register(MulticoreParam(workers=n_cores))
    registerDoParallel(cores=n_cores, )
    message("Computing person_variant.")

    # Compute p-values.
    message("Computing p-values")
    if (paired) {
        permuted_pairs = PermuteDataset(ids = group_idx[which(phenotypes)], n_permute_iter = n_permutations)
        p_all = data.frame(
            foreach(i_group=1:ncol(person_variants))
            %dopar% ComputePairedP(person_variants, phenotypes = phenotypes, shuff_switch_sign = permuted_pairs, i_group = i_group))
        }
    else {
        permuted_phenos = PermuteDataset(ids = group_idx, n_permute_iter = n_permutations, prob=c(mean(phenotypes==T), mean(phenotypes==F)))
        p_all = data.frame(
            foreach(i_group=1:ncol(person_variants))
            %dopar% ComputeUnpairedP(person_variants, phenotypes = phenotypes, shuff_phenos = permuted_phenos, i_group = i_group))
    }
    
    # Set up exports.
    p = melt(p_all[1,])$value # P-values for the true data.
    p_shuff = p_all[-1,] # P-values for our shuffled data.
    t = SetTime(t=NA)
    message("Exporting Data.")
    export = list()
    export$p = p
    export$p_shuff = p_shuff
    export$var_names = var_names
    export$person_variants = person_variants
    if (paired) { export$permuted_pairs = permuted_pairs }
    else {export$permuted_phenos = permuted_phenos }
    return(export)
    message("Success.")
}
