library(Matrix)
source("/oak/stanford/groups/dpwall/users/briannac/phyloWAS/R/helpers.R")

#' Turns matrix of taxa vs variants into a sparse matrix.
#'
#' @param taxa_vars matrix of taxa vs variants.
#' @param variant_vals list of variant values to check for.
#' @return sparse representation of taxa variant matrix.
#'
#' @export
TaxaVariantSparse = function(taxa_vars, variant_vals) {
    var_names = c()
    if (length(variant_vals)==0) {
        taxa_vars_sparse = as.sparseMatrix(taxa_vars)
        var_names = 1:length(taxa_vars)
    } else {
        taxa_vars_sparse = sparseMatrix(c(1), c(1),
                                        dims = c(nrow(taxa_vars), ncol(taxa_vars)*length(variant_vals)))
        for (i_variant in 1:length(variant_vals)) {
            v = variant_vals[i_variant]
            ids = (ncol(taxa_vars)*(i_variant-1)+1):(ncol(taxa_vars)*(i_variant))
            taxa_vars_sparse[,ids] = as.sparseMatrix(taxa_vars==v)
            var_names = c(var_names, sapply(1:ncol(taxa_vars), function(pos) {paste0(pos, '.', v)})) 
        }
    }
    row.names(taxa_vars_sparse) = row.names(taxa_vars)
    colnames(taxa_vars_sparse) = var_names

    return(taxa_vars_sparse)
}


#' Computes combinations of n_df variants.
#'
#' @param taxa_variant_sparse sparse matrix of taxa vs variants.
#' @param i_taxa which taxa to compute variants for.
#' @param n_df number of variants (combos of 2, 3, 4, etc.) to consider in each combination.
#' @param index_lookup an n_df matrix containing the element index for each n_df index (row, col, ...) 
#' @return sparse representation of taxa variant matrix.
#'
#' @export vector containing info for which variant combinations occur in i_taxa.
TaxaVariantCombos = function(taxa_variant_sparse, i_taxa, n_df) {
    
    # Find which variants actually exist in a given taxa.
    # Create a sparse matrix with n_df columns and n_variants rows
    # Each column will have same value (1 or 0 if variant exists.)
    var_ids=which(taxa_variant_sparse[i_taxa,]!=0)
    for (i in 1:n_df) {
        if (i==1) {df = data.frame(var_ids)}
        else {df = cbind(df, var_ids)}
    }
    # expand.grid finds all the combinations of variants 
    # sort of does 
    # (all possible combos of (x,y) for x in col1, and y in col2).
    idx = expand.grid(df)
    if ((i_taxa %% 100) == 0) {
        message(paste0("On taxa ", i_taxa))
    }
    #colnames(idx) = 1:n_df
    if (n_df >= 2) {
        for (i in 2:n_df) {
            idx = idx[idx[,i-1] < idx[,i],]
        }
    }

    idx = as.matrix(idx)
    
    #element_ids = lookup_table[idx]
    element_ids = ElementIdsFromDimensionIds(idx, rep(ncol(taxa_variant_sparse), n_df))
    mat = sparseMatrix(i = rep(1, length(element_ids)), j = element_ids,dims = c(1,ncol(taxa_variant_sparse)**n_df))
    return(mat)
}
   
                          
#' Computes the abundance of each 'variant' for each person
#'
#' @param person_taxa matrix of person vs taxa abundance.
#' @param taxa_var matrix of taxa vs variant (presence T/F)
#' @return sparse representation of abundance of each variant for each person.
#'
#' @export
PersonVars = function(person_taxa, taxa_var) {
    return(person_taxa %*% taxa_var)
}


#' Computes the abundance of each 'variant' for each person
#'
#' This function trims the a matrix to only include non-duplicated columns.
#'
#' @param taxa_var matrix of person vs taxa abundance.
#' @param has_vals whether or not the sparse representation contains values (T) or if it just contains 0/1s (F).
#' @return trimmed 
#'
#' @export
TrimVarMatrix = function(mat, has_vals=F) {
    # Get rid of variants that don't exist in any taxa.
    #t=Sys.time()
    #to_keep = which(apply(mat>0, 2, any))
    #mat = mat[,to_keep]
    #var_names = var_names[to_keep]
    #t = SetTime(t)
    message("trimvarmatrix 1")
    # Get rid of variants that are perfectly in LD.
    to_keep = which(!duplicated.sparseMatrix(mat, 2, has_vals=has_vals))
    #mat = mat[,to_keep]
    #var_names = var_names[to_keep]
    #message("trimvarmatrix 2")
    #t = SetTime(t)

    # Get rid of variants that only exist in a single taxa.
    #to_keep = unique(Matrix::summary(taxa_var_trimmed)[,'j'])
    #taxa_var_trimmed = taxa_var_trimmed[, to_keep]
    #var_names_trimmed = var_names_trimmed[to_keep]
    return(to_keep)
}