context("helpers")

test_that("GetPrimes works", {
    expect_equal(GetPrimes(5), c(2,3,5))
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


test_that("duplicated.sparseMatrix works", {
    mat =  sparseMatrix(c(1,1,2, 2,2,3,3), c(1,2,1, 2,3,1,2), x=c(1,1,1,1,1,1,1), dims=c(3,3))
    expect_equal(duplicated.sparseMatrix(mat, 1), c(F,F,T))
    expect_equal(duplicated.sparseMatrix(mat, 2), c(F,T,F))
    expect_equal(duplicated.sparseMatrix(mat, 1, has_vals=T), c(F,F,T))
    expect_equal(duplicated.sparseMatrix(mat, 2, has_vals=T), c(F,T,F))
})

