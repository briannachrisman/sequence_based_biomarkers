context("statistics")

test_that("WilcoxonRankFastPaired works", {
    #expect_equal(WilcoxonRankFastPaired(c(1,2,5,3,4), c(F,F,T,T,T)), 0.15625)
    phenotypes = c(T,F,F,T,F,T,T,F,T,F)
    values = c(.2,.4,.5,.6,.6, .6,.3,.4,.3,.2)
    wilcox_r = wilcox.test(values[phenotypes==T], values[phenotypes==F], paired = T, exact=F, alternative='two.sided')$p.value
    diffs = values[phenotypes==T] - values[phenotypes==F]
    diffs = diffs[diffs!=0]
    signs = diffs>0
    ranks = rank(abs(diffs))
    expect_equal(WilcoxonRankFastPaired(ranks, signs), wilcox_r)
})

test_that("WilcoxonRankFastUnpaired works", {
    phenotypes = c(T,T,T,T,F,F,F,F,T,T)
    values = c(.2,.4,.5,.6,.6,.6,.3,.4,.3,.2)
    wilcox_r = wilcox.test(values[phenotypes==T], values[phenotypes==F], paired = F, exact=F, alternative='two.sided')$p.value
    ranks = rank(values)
    wilcox_fast = WilcoxonRankFastUnpaired(ranks, phenotypes)
    expect_equal(wilcox_r, wilcox_fast)
})

test_that("ComputePairedP works", {
    person_groups = data.frame(group1=c(.01,.03,0,0,.06,.1, .2, .4, .3, .9, .31, .99), group2=c(.1,.2,0,0,0,0,.2,.1,.1,.2,.2,.1))
    phenotypes = c(T,F,T,F,T,F,T,F,T,F,T,F)
    shuff_switch_sign = data.frame(c(T,T,T,T,T,T), c(T,F,F,F,T,T))
    p = ComputePairedP(person_groups, phenotypes = phenotypes, shuff_switch_sign, 1)
    expect_equal(round(p, 3), c(0.059, 0.059, 0.590))
})

test_that("ComputeUnpairedP works", {
    person_groups = data.frame(group1=c(.01,.02,.05,.1,.2,.5, .6, .2, .8, .9, .88, .88), group2=c(.1,.2,0,0,0,0,.2,.1,.1,.2,.2,.1))
    phenotypes = c(T,F,F,T,T,F,F,F,T,F,F,F)
    shuff_switch_sign = data.frame(c(T,F,F,T,T,F,F,F,T,F,F,F), c(T,F,T,F,F,T,T,T,F,F,F,T))
    p = ComputeUnpairedP(person_groups, phenotypes = phenotypes, shuff_switch_sign, i_group=1)
    expect_equal(round(p, 3), c(0.306, 0.306, 0.574))
})

test_that("PermuteDataset works", {
    set.seed(42)
    expect_equal(PermuteDataset(c(1,2,2,3,3,4,7), 2),matrix(c(c(T,T,T,F,F,T,T), c(T,T,T,F,F,T,T)), ncol=2))
})