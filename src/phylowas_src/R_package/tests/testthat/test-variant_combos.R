context("variant_combos")

test_that("TaxaVariantSparse works", {
    taxa_vars = t(data.frame(ASV1=c('A', 'T', 'G', '-', '-', 'X'),
                         ASV2=c('T', '-', 'T', '-', 'C', 'X')))
    variant_vals = c('A','T','G','C')
    variant_sparse_mat = TaxaVariantSparse(taxa_vars, variant_vals)
    expected = sparseMatrix(c(1,1,1,2,2,2), c(1,8, 15, 7, 9, 23), dims=c(2,24))
    rownames(expected) = c('ASV1', 'ASV2')
    colnames(expected) = c('1.A', '2.A', '3.A', '4.A', '5.A', '6.A',
                           '1.T', '2.T', '3.T', '4.T', '5.T', '6.T',
                           '1.G', '2.G', '3.G', '4.G', '5.G', '6.G',
                           '1.C', '2.C', '3.C', '4.C', '5.C', '6.C')
    expect_equal(variant_sparse_mat, expected)
})


test_that("TaxaVariantCombos works", {
    taxa_variant_sparse = sparseMatrix(c(1,1,2,2), c(1,2, 4, 4), dims=c(2,4))
    n_df = 2
    idx_lookup_table =  array(1:(ncol(taxa_variant_sparse)**n_df), rep(ncol(taxa_variant_sparse), n_df))
    expect_equal(TaxaVariantCombos(taxa_variant_sparse, 1, 2, idx_lookup_table), data.frame(j=c(1,2,5,6),i=1))
    expect_equal(TaxaVariantCombos(taxa_variant_sparse, 2, 2, idx_lookup_table), data.frame(j=c(16),i=2))
})


test_that("TrimTaxaVar works", {
    mat =  sparseMatrix(c(1,1,2, 2,2,3,3), c(1,2,1, 2,3,1,2), x=c(1,1,1,1,1,1,1), dims=c(3,4))
    var_names = c('var1', 'var2', 'var3', 'var4')
    trimmed=TrimTaxaVar(mat, var_names)$taxa_var_trimmed
    trimmed_var_names=TrimTaxaVar(mat, var_names)$var_names_trimmed
    expected = sparseMatrix(c(1,2,2,3), c(1,1,2,1), x=1, dims=c(3,2))
    expect_equal(trimmed, expected)
    expect_equal(trimmed_var_names, c('var1', 'var3'))
})