# This tests the correctGroupSummary function.
# library(testthat); library(scuttle); source("test-correct-summary.R")

library(scuttle)

test_that("correctGroupSummary works correctly in raw mode", {
    # No batches.
    y <- matrix(rnorm(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged <- correctGroupSummary(y, group, block)
    ref <- sumCountsAcrossCells(y, group, average=TRUE)
    expect_equal(averaged, assay(ref))

    averaged <- correctGroupSummary(y, group, block=NULL) # same with NULL.
    expect_equal(averaged, assay(ref))
 
    # Perfectly balanced.
    y <- matrix(rnorm(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep(1:10, 10)
    averaged <- correctGroupSummary(y, group, block)
 
    ref <- sumCountsAcrossCells(y, group, average=TRUE)
    ref <- assay(ref)
 
    ref2 <- ref - rowMeans(ref)
    averaged2 <- averaged - rowMeans(averaged)
    expect_equal(ref2, averaged2)
    
    # Effectively ignores batches with only one cluster.
    y <- matrix(rnorm(210), ncol=21)
    group <- c(1:10, 1:10, 1)
    block <- rep(1:3, c(10, 10, 1))
    averaged <- correctGroupSummary(y, group, block)
    averaged2 <- correctGroupSummary(y[,-21], group[-21], block[-21])
    expect_equal(averaged, averaged2)
 
    # Handles batch-specific clusters.
    y <- matrix(rnorm(220), ncol=22)
    group <- c(1:10, 1:10, c(1,11))
    block <- rep(1:3, c(10, 10, 2))
    averaged <- correctGroupSummary(y, group, block)
    expect_identical(sum(!colAnyNAs(averaged)), 11L)
})

test_that("correctGroupSummary works with weights", {
    y <- matrix(rnorm(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep(1:10, 10)

    w <- sample(2, ncol(y), replace=TRUE)
    out <- correctGroupSummary(y, group, block, weights=w)

    expander <- rep(seq_along(w), w)
    ref <- correctGroupSummary(y[,expander], group[expander], block[expander])
    expect_equal(out, ref)
})

test_that("correctGroupSummary handles subsetting", {
    y <- matrix(rnorm(1000), ncol=100)
    group <- factor(rep(1:10, each=10))
    block <- rep(1:10, 10)

    # Subset.row gives the same results.
    sub <- sample(nrow(y), 10)
    ref <- correctGroupSummary(y, group, block, subset.row=sub)
    out <- correctGroupSummary(y[sub,,drop=FALSE], group, block)
    expect_identical(ref, out)

    # Handles the special one-gene case.
    ref <- correctGroupSummary(y, group, block)
    out <- correctGroupSummary(y[1,,drop=FALSE], group, block)
    expect_identical(ref[1,,drop=FALSE], out)

    w <- runif(100)
    ref <- correctGroupSummary(y, group, block, weights=w)
    out <- correctGroupSummary(y[1,,drop=FALSE], group, block, weights=w)
    expect_identical(ref[1,,drop=FALSE], out)
})

test_that("correctGroupSummary respects factor ordering", {
    y <- matrix(rnorm(1000), ncol=100)
    group <- factor(rep(1:10, each=10), 10:1)
    block <- rep(1:10, 10)

    out <- correctGroupSummary(y, group, block)
    expect_identical(colnames(out), as.character(10:1))

    ref <- correctGroupSummary(y, as.character(group), block)
    expect_equal(out, ref[,colnames(out)])
})

test_that("correctGroupSummary works correctly in log mode", {
    # Actually has an effect.
    y <- matrix(rexp(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged1 <- correctGroupSummary(y, group, block)
    averaged2 <- correctGroupSummary(y, group, block, transform="log")
    expect_false(identical(averaged1, averaged2))

    # Survives a round trip with correct untransformation. 
    y <- matrix(rexp(1000), ncol=10)
    group <- 1:10
    block <- rep('all', ncol(y))
    averaged <- correctGroupSummary(y, group, block, transform="log")
    expect_equivalent(averaged, y)

    # Handles zeroes.
    y <- matrix(0, ncol=100, nrow=10)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged <- correctGroupSummary(y, group, block, transform="log")
    expect_true(all(averaged==0))
})

test_that("correctGroupSummary works correctly in logit mode", {
    # Actually has an effect.
    y <- matrix(runif(1000), ncol=100)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged1 <- correctGroupSummary(y, group, block)
    averaged2 <- correctGroupSummary(y, group, block, transform="log")
    averaged3 <- correctGroupSummary(y, group, block, transform="logit")
    expect_false(identical(averaged1, averaged3))
    expect_false(identical(averaged2, averaged3))

    # Survives a round trip with correct untransformation. 
    y <- matrix(runif(1000), ncol=10)
    group <- 1:10
    block <- rep('all', ncol(y))
    averaged <- correctGroupSummary(y, group, block, transform="logit")
    expect_equivalent(averaged, y)

    # Handles boundaries.
    y <- matrix(0, ncol=100, nrow=10)
    group <- rep(1:10, each=10)
    block <- rep('all', ncol(y))
    averaged <- correctGroupSummary(y, group, block, transform="logit")
    expect_true(all(abs(averaged) < 1e-10))

    y <- matrix(1, ncol=100, nrow=10)
    averaged <- correctGroupSummary(y, group, block, transform="logit")
    expect_true(all(averaged==1))
})

