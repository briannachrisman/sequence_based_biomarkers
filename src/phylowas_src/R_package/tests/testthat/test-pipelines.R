context("pipelines")

test_that("RunSingleVariantMWAS works", {
person_taxa = sparseMatrix(c(1,1,2,2,3,3,4,4,5,5,6,6), c(1,2,1,3,1,4,2,5,5,6,5,6), x=c(.1,.3,.2,4,.6,.5,.4,.5,.5,.5,.3,.3))
    taxa_variant_combos = t(data.frame(ASV1=c(T,T,T,F,F),
                        ASV2=c(T,T,T,F,F),
                        ASV3=c(T,F,T,T,F),
                        ASV4=c(F,F,F,F,T),
                        ASV5=c(T,T,F,F,F),
                        ASV6=c(T,T,T,T,T)))
    phenotypes = c(T,F,T,F,T,F)
    idx = c(1,1,2,2,4,4)
    set.seed(43)
    expect_equal(is.list(RunSingleVariantMWAS(person_taxa, taxa_variant_combos, phenotypes, idx, n_df, paired=F, n_permutations=2, n_cores=1)),
                 T)
    expect_equal(is.list(RunSingleVariantMWAS(person_taxa, taxa_variant_combos, phenotypes, idx, n_df, paired=T, n_permutations=2, n_cores=1)),
                 T)
})


test_that("RunVariantComboMWAS works", {
    person_taxa = sparseMatrix(c(1,1,2,2,3,3,4,4,5,5,6,6), c(1,2,1,3,1,4,2,5,5,6,5,6), x=c(.1,.3,.2,4,.6,.5,.4,.5,.5,.5,.3,.3))
    seq_df = t(data.frame(ASV1=c('A','X','T','-','X'),
                        ASV2=c('T','X','C','-','X'),
                        ASV3=c('T','G','T','-','X'),
                        ASV4=c('T','X','C','T','X'),
                        ASV5=c('X','G','T','-','X'),
                        ASV6=c('-','C','C','T','X')))
    phenotypes = c(T,F,T,F,T,F)
    idx = c(1,1,2,2,4,4)
    paired=T
    bases = c('A','C','T', 'G')
    set.seed(43)
    ans_unpaired = RunVariantComboMWAS(person_taxa = person_taxa, seq_df = seq_df, phenotypes = phenotypes, idx = idx, n_df = 2, paired = T, bases = bases, n_permutations=2, n_cores=1)
    expect_equal(is.list(ans_unpaired), T)
    
    # Check if tables line up.
    ps = ans_unpaired$ps
    p = ans_unpaired$p
    interesting_taxa = taxa_names(ans_unpaired$ps)
    variants = strsplit(import$var_names[order(p)[i]], split = '_')[[1]]
    for (v in variants) {
        print(v)
        var_split = strsplit(v, split = '[.]')
        pos = as.integer(var_split[[1]][1])
        base = var_split[[1]][2]
        interesting_taxa = intersect(interesting_taxa, names(which(ans_unpaired$seq_df[,pos]==base)))
    }
    expect_equal(sum(otu_table(ps)[,interesting_taxa]), sum(ans_unpaired$person_variants[, order(p)[i]]))
    
    
})

test_that("SetTime works", {
    t1  = Sys.time()
    Sys.sleep(1)
    t2 = SetTime(t1)
    expect_equal(round(as.numeric(t2-t1)), 1)
})

test_that("as.sparseMatrix works", {
    mat = matrix(c(c(0,2,0), c(1,2,0)), ncol = 2)
    expect_equal(as.sparseMatrix(mat, use_pointers = T),
                 sparseMatrix(c(1,2,2), c(2,1,2), x=c(1,2,2), dims = dim(mat)))
    expect_equal(as.sparseMatrix(mat),
                 sparseMatrix(c(1,2,2), c(2,1,2), dims = dim(mat)))
})
