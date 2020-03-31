#' Computes the p-value of the wilcoxon rank test.
#'
#' This function performs a wilcox rank sum test given the ranking and the signs.
#'
#' @param ranks list of rank differences.
#' @param signs list sof sign differences.
#' @return p-value.
#'
#' @export
WilcoxonRankFastPaired = function(ranks, signs) {
    #n = length(ranks)
    #W_T = abs(sum(ranks[signs==T]))
    #W_F = abs(sum(ranks[signs==F]))
    #return(psignrank(min(W_T, W_F), n=n))

    n = length(ranks)
    n_ties <- table(ranks)
    W = sum(ranks[signs > 0])
    z = W - n * (n + 1)/4
    sigma = sqrt(n * (n + 1) * (2 * n + 1) / 24
                              - sum(n_ties^3 - n_ties) / 48)
    correction = sign(z) * 0.5
    z = (z - correction) / sigma
    return(2 * min(pnorm(z), pnorm(z, lower.tail=FALSE)))

}

WilcoxonRankFastUnpaired = function(ranks, signs) {
    x = ranks[signs==T]
    y = ranks[signs==F]
    r = c(x, y)
    n.x = length(x)
    n.y = length(y)
    W  = sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2
    n_ties = table(r)
    z = W - n.x * n.y / 2
    sigma = sqrt((n.x * n.y / 12) *
                              ((n.x + n.y + 1)
                               - sum(n_ties^3 - n_ties)
                               / ((n.x + n.y) * (n.x + n.y - 1))))
    correction = sign(z) * 0.5
    z = (z - correction) / sigma
    return(2 * min(pnorm(z), pnorm(z, lower.tail=FALSE)))
}

#' Computes the p-value for a paired wilcoxon rank test.
#'
#' This function computes the p-value for a paired wilcoxon rank sum test, and computes p-values from the permuted data.
#'
#' @param person_groups table of person (rows) vs abundance (columns) for groups of taxa or variants.
#' @param phenotypes list of phenotypes for each person (T/F).
#' @param shuff_switch_sign table of pair (row) vs iteration (col) of whether or not to switch the phenotypes of each pair for that permutation.
#' @param i_group which group (col) of person_groups you would like to compute the p-value for.
#' @return p-value.
#'
#' @export
ComputePairedP = function(person_groups, phenotypes, shuff_switch_sign, i_group)  {
    ids_T = which(phenotypes==T) # Rows ids for which the phenotype==T.
    ids_F = which(phenotypes==F) # Row ids for which phenotype==F.
    n_shuff_iter = ncol(shuff_switch_sign)
    # Only look at pairs in which neither person has an abundance of 0 for group i_group.
    ids = union(which(person_groups[ids_T,i_group]!=0), which(person_groups[ids_F,i_group]!=0))
        
    # Compute the difference in abundance between corresponding pairs, transform to ranks, compute signs.
    diffs = person_groups[ids_T,i_group][ids] - person_groups[ids_F,i_group][ids]
    shuff_switch_sign = shuff_switch_sign[ids, ]
    if (sum(diffs!=0)<2) {return(append(1, sapply(1:n_shuff_iter, function(x) {1})))} # If not enough diffs, then p-value = 1.
    diffs_ranks = rank(abs(diffs))
    diffs_signs = ifelse(diffs>0, yes=T, no=F)
                                                      

    # Compute p given the current data.
    p = WilcoxonRankFastPaired(diffs_ranks, diffs_signs)
    # Compute p-value for every shuffled set of phenotypes for i_group.
    p_shuffs = sapply(1:n_shuff_iter,
                      function(i_iter) {WilcoxonRankFastPaired(diffs_ranks, diffs_signs!=shuff_switch_sign[,i_iter])})  # Compute p-values of shuffled datasets
  
    # Print progress every ~1% increase.
    if ((i_group/ncol(person_groups)*100) %% 10 < .01) {
        print(paste0(100*i_group/ncol(person_groups), "% done with p-value computation."))
        flush.console()
    }
    return(append(p, p_shuffs))
}

#' Computes the p-value for an unpaired wilcoxon rank test.
#'
#' This function computes the p-value for an unpaired wilcoxon rank sum test, and computes p-values from the permuted data.
#'
#' @param person_groups table of person (rows) vs abundance (columns) for groups of taxa or variants.
#' @param phenotypes list of phenotypes for each person (T/F).
#' @param shuff_switch_sign table of person (row) vs iteration (col) of whether or not to switch the phenotypes of each person for that permutation.
#' @param i_group which group (col) of person_groups you would like to compute the p-value for.
#' @return p-value.
#'
#' @export
ComputeUnpairedP = function(person_groups, phenotypes, shuff_phenos, i_group)  {
    n_shuff_iter = ncol(shuff_phenos)
    
    # Compute rankings and signs to pass to wilcoxon test.
    vals = person_groups[,i_group]
    if (sum(vals!=0)<2) {return(append(1, sapply(1:n_shuff_iter, function(x) {1})))} # If not enough diffs, then p-value = 1.
    vals_ranks = rank(vals)
    vals_signs = phenotypes
                                                      
    # Compute p given the current data.
    p = WilcoxonRankFastUnpaired(vals_ranks, vals_signs)
    
    # Compute p-value for every shuffled set of phenotypes for i_group.
    p_shuffs = sapply(1:n_shuff_iter,
                      function(i_iter) {WilcoxonRankFastUnpaired(vals_ranks, sample(phenotypes))})  # Compute p-values of shuffled datasets
  
    # Print progress every ~1% increase.
    if ((i_group/ncol(person_groups)*100) %% 10 < .01) {
        print(paste0(100*i_group/ncol(person_groups), "% done with p-value computation."))
    }
    return(append(p, p_shuffs))
}

#' Computes a matrix illustrating how to permute the dataset.
#'
#' This function creates a matrix with length(ids) rows and n_permute_iter columns, showing how to permute the dataset
#' for each permutation (whether or not to flip the phenotype at row (row # for permutation #).
#'
#' @param ids identifier of rows you want to keep together. IE if you have a family, and always want to flip their phenotypes together, those two rows should have the same ID.
#' @param n_permute_iter number of permutations you want to run.
#' @return permute_switch_pheno a matrix with length(ids) rows and n_permute_iter columns, showing how to permute the dataset for each permutation.TrimTaxaVar
#'
#' @export
PermuteDataset = function(ids, n_permute_iter, prob=c(.5,.5))  {
    factor_ids = factor(ids)
    permute_switch_pheno = array(dim = c(length(factor_ids), n_permute_iter))
    for (i_permute in 1:n_permute_iter) {
        id_table = sapply(unique(factor_ids), function(x){sample(c(T, F), size=1, prob=prob)})
        names(id_table) = unique(factor_ids)
        permute_switch_pheno[,i_permute] = sapply(factor_ids, function(x) {id_table[x]})
    }
    return (permute_switch_pheno)
}


