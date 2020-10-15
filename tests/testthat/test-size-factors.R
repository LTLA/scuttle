# Tests for size factor calculations.
# library(scuttle); library(testthat); source("setup.R"); source("test-size-factors.R")

X <- sce
dummy <- counts(X)

test_that("librarySizeFactors works as expected", {
    sf <- librarySizeFactors(X)
    expect_identical(mean(sf), 1)
    expect_true(sd(sf/colSums(dummy)) < 1e-8)

    X <- computeLibraryFactors(X)
    expect_identical(sf, sizeFactors(X))

    sf <- librarySizeFactors(X, subset.row=1:10)
    expect_identical(mean(sf), 1)
    expect_true(sd(sf/colSums(dummy[1:10,])) < 1e-8)
})

test_that("geometric size factors work as expected", {
    sf <- geometricSizeFactors(X)
    expect_equal(mean(sf), 1)
    expect_true(sd(sf/exp(colMeans(log1p(dummy)))) < 1e-8)

    X <- computeGeometricFactors(X)
    expect_identical(sf, sizeFactors(X))

    sf <- geometricSizeFactors(X, subset.row=100:200)
    expect_identical(mean(sf), 1)
    expect_true(sd(sf/exp(colMeans(log1p(dummy[100:200,])))) < 1e-8)
})

test_that("medianSizeFactors works as expected", {
    sf <- medianSizeFactors(X)
    expect_equal(mean(sf), 1)

    X <- computeMedianFactors(X)
    expect_identical(sf, sizeFactors(X))

    sf <- medianSizeFactors(X, subset.row=1:10)
    expect_equal(mean(sf), 1)
    expect_identical(sf, medianSizeFactors(X[1:10,]))

    ref <- runif(nrow(X))
    sf <- medianSizeFactors(X, subset.row=100:200, reference=ref)
    expect_identical(mean(sf), 1)
    expect_identical(sf, medianSizeFactors(X[100:200,], reference=ref[100:200]))
})
