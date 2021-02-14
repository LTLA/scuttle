# Testing the downsampling functions.
# library(scuttle); library(testthat); source("test-downsample.R")

CHECKFUN <- function(input, prop) {
    out <- downsampleMatrix(input, prop)
    expect_identical(colSums(out), round(colSums(input)*prop))
    expect_true(all(out <= input))
    return(invisible(NULL))
}

CHECKSUM <- function(input, prop) {
    out <- downsampleMatrix(input, prop, bycol=FALSE) 
    expect_equal(sum(out), round(prop*sum(input)))
    expect_true(all(out <= input))
    return(invisible(NULL))
}

test_that("downsampling from a count matrix gives expected sums", {
    # Vanilla run.
    set.seed(0)
    ncells <- 100
    u1 <- matrix(rpois(20000, 5), ncol=ncells)
    u2 <- matrix(rpois(20000, 1), ncol=ncells)

    set.seed(100)
    for (down in c(0.111, 0.333, 0.777)) { # Avoid problems with different rounding of 0.5.
        CHECKFUN(u1, down) 
        CHECKSUM(u1, down) 
    }

    set.seed(101)
    for (down in c(0.111, 0.333, 0.777)) { # Avoid problems with different rounding of 0.5.
        CHECKFUN(u2, down) 
        CHECKSUM(u2, down) 
    }

    # Checking double-precision inputs.
    v1 <- u1
    storage.mode(v1) <- "double"
    set.seed(200)
    for (down in c(0.111, 0.333, 0.777)) { 
        CHECKFUN(v1, down) 
        CHECKSUM(v1, down) 
    }

    v2 <- u2
    storage.mode(v2) <- "double"
    set.seed(202)
    for (down in c(0.111, 0.333, 0.777)) { 
        CHECKFUN(v2, down) 
        CHECKSUM(v2, down) 
    }

    # Checking vectors of proportions.
    set.seed(300)
    CHECKFUN(u1, runif(ncells))
    CHECKFUN(u1, runif(ncells, 0, 0.5))
    CHECKFUN(u1, runif(ncells, 0.1, 0.2))

    set.seed(303)
    CHECKFUN(u2, runif(ncells))
    CHECKFUN(u2, runif(ncells, 0, 0.5))
    CHECKFUN(u2, runif(ncells, 0.1, 0.2))

    # Checking that bycol=FALSE behaves consistently with bycol=TRUE. 
    set.seed(505)
    out1 <- downsampleMatrix(u1, prop=0.111, bycol=FALSE)
    set.seed(505)
    ref <- downsampleMatrix(cbind(as.vector(u1)), prop=0.111, bycol=TRUE)
    dim(ref) <- dim(out1)
    expect_identical(ref, out1)
})

test_that("downsampling from a count matrix worsk with silly inputs", {
    ncells <- 100
    u1 <- matrix(rpois(20000, 5), ncol=ncells)
    expect_equivalent(as.matrix(downsampleMatrix(u1[0,,drop=FALSE], prop=0.5)), u1[0,,drop=FALSE])
    expect_equivalent(as.matrix(downsampleMatrix(u1[,0,drop=FALSE], prop=0.5)), u1[,0,drop=FALSE])
    expect_equivalent(as.matrix(downsampleMatrix(u1[0,0,drop=FALSE], prop=0.5)), u1[0,0,drop=FALSE])

    v1 <- u1
    storage.mode(v1) <- "double"
    expect_equivalent(as.matrix(downsampleMatrix(u1[0,,drop=FALSE], bycol=TRUE, prop=0.5)), u1[0,,drop=FALSE])
    expect_equivalent(as.matrix(downsampleMatrix(u1[,0,drop=FALSE], bycol=TRUE, prop=0.5)), u1[,0,drop=FALSE])
    expect_equivalent(as.matrix(downsampleMatrix(u1[0,0,drop=FALSE], bycol=TRUE, prop=0.5)), u1[0,0,drop=FALSE])

    w1 <- as(u1, "dgCMatrix")
    expect_equivalent(downsampleMatrix(w1[0,,drop=FALSE], prop=0.5), w1[0,,drop=FALSE])
    expect_equivalent(downsampleMatrix(w1[,0,drop=FALSE], prop=0.5), w1[,0,drop=FALSE])
    expect_equivalent(downsampleMatrix(w1[0,0,drop=FALSE], prop=0.5), w1[0,0,drop=FALSE])
})

test_that("different matrix representations yield the same result", {
    set.seed(500)
    ncells <- 100
    u1 <- matrix(rpois(20000, 5), ncol=ncells)
    v1 <- as(u1, "dgCMatrix")
    w1 <- as(u1, "dgTMatrix")

    # Basic downsampling.
    for (down in c(0.111, 0.333, 0.777)) { 
        set.seed(501)
        dd <- downsampleMatrix(u1, down)
        expect_s4_class(dd, "dgCMatrix")

        set.seed(501)
        dc <- downsampleMatrix(v1, down)
        expect_identical(dc, dd)

        set.seed(501)
        dt <- downsampleMatrix(w1, down)
        expect_identical(dt, dd)
    }

    # Columnar downsampling.
    for (down in c(0.111, 0.333, 0.777)) { 
        set.seed(502)
        dd <- downsampleMatrix(u1, down, bycol=TRUE)
        expect_s4_class(dd, "dgCMatrix")

        set.seed(502)
        dc <- downsampleMatrix(v1, down, bycol=TRUE)
        expect_identical(dc, dd)

        set.seed(502)
        dt <- downsampleMatrix(w1, down, bycol=TRUE)
        expect_identical(dt, dd)
    }

    # Columnar downsampling.
    prop <- runif(ncol(u1))

    set.seed(503)
    dd <- downsampleMatrix(u1, prop, bycol=TRUE)
    expect_s4_class(dd, "dgCMatrix")

    set.seed(503)
    dc <- downsampleMatrix(v1, prop, bycol=TRUE)
    expect_equivalent(dc, dd)

    set.seed(503)
    dt <- downsampleMatrix(w1, prop, bycol=TRUE)
    expect_equivalent(dt, dd)
})

set.seed(510)
test_that("downsampleMatrix responds to various DelayedArray options", {
    ncells <- 100
    u1 <- matrix(rpois(20000, 5), ncol=ncells)
    prop <- runif(ncol(u1))

    set.seed(504)
    refF <- downsampleMatrix(u1, 0.211, bycol=FALSE)
    set.seed(504)
    refT <- downsampleMatrix(u1, prop, bycol=TRUE)

    library(DelayedArray)
    D1 <- DelayedArray(u1)
    old <- getAutoBlockSize()
    for (block.size in c(1000, 10000, 100000)) {
        setAutoBlockSize(block.size)

        set.seed(504)
        obsF <- downsampleMatrix(D1, 0.211, bycol=FALSE)
        expect_identical(refF, obsF)

        set.seed(504)
        obsT <- downsampleMatrix(D1, prop, bycol=TRUE)
        expect_identical(refT, obsT)
    }

    setAutoBlockSize(old)

    # Setting the realization sink.
    sink <- AutoRealizationSink(dim(u1)) 
    set.seed(504)
    sunkF <- downsampleMatrix(u1, 0.211, bycol=FALSE, sink=sink)
    expect_s4_class(sunkF, "DelayedMatrix")
    expect_identical(unname(as.matrix(refF)), as.matrix(sunkF))

    sink <- AutoRealizationSink(dim(u1)) 
    set.seed(504)
    sunkT <- downsampleMatrix(u1, prop, bycol=TRUE, sink=sink)
    expect_s4_class(sunkT, "DelayedMatrix")
    expect_identical(unname(as.matrix(refT)), as.matrix(sunkT))
})

set.seed(500)
test_that("downsampling from a count matrix gives expected margins", {
    # Checking that the sampling scheme is correct (as much as possible).
    known <- matrix(1:5, nrow=5, ncol=10000)
    prop <- 0.51
    truth <- known[,1]*prop
    out <- downsampleMatrix(known, prop)
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.1)) # Less than 10% error on the estimated proportions.

    out <- downsampleMatrix(known, prop, bycol=FALSE) # Repeating by column.
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.1)) 

    # Repeating with larger counts.
    known <- matrix(1:5*100, nrow=5, ncol=10000)
    prop <- 0.51
    truth <- known[,1]*prop
    out <- downsampleMatrix(known, prop)
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.01)) # Less than 1% error on the estimated proportions.

    out <- downsampleMatrix(known, prop, bycol=FALSE)
    expect_true(all(abs(rowMeans(out)/truth - 1) < 0.01)) 

    # Checking the column sums when bycol=FALSE.
    known <- matrix(100, nrow=1000, ncol=10)
    out <- downsampleMatrix(known, prop, bycol=FALSE)
    expect_true(all(abs(colMeans(out)/colMeans(known)/prop - 1) < 0.01))

    # Checking that downsampling preserves relative abundances.
    set.seed(500)
    X <- matrix(1:4*100, ncol=500, nrow=4)
    Y <- downsampleMatrix(X, prop=0.11)
    expect_true(all(abs(rowMeans(Y) - 0.11*rowMeans(X)) < 1))
    Y <- downsampleMatrix(X, prop=0.55)
    expect_true(all(abs(rowMeans(Y) - 0.55*rowMeans(X)) < 1))
    Y <- downsampleMatrix(X, prop=0.11, bycol=FALSE)
    expect_true(all(abs(rowMeans(Y) - 0.11*rowMeans(X)) < 1))
    Y <- downsampleMatrix(X, prop=0.55, bycol=FALSE)
    expect_true(all(abs(rowMeans(Y) - 0.55*rowMeans(X)) < 1))
})

set.seed(5001)
test_that("downsampling batches gives expected results", {
    u1 <- matrix(rpois(20000, 5), ncol=100)

    set.seed(0)
    output <- downsampleBatches(u1, u1*10)
    expect_equal(colSums(output[[1]])/colSums(output[[2]]), rep(1, ncol(output[[1]])))
    expect_equivalent(u1, as.matrix(output[[1]]))

    # Checking that the output responds to the seed.
    set.seed(0)
    output2 <- downsampleBatches(u1, u1*10)
    expect_identical(output, output2)
    output2 <- downsampleBatches(u1, u1*10)
    expect_false(identical(output, output2))

    # Works with bycol=FALSE.
    output <- downsampleBatches(u1, u1*10, bycol=FALSE)
    expect_equal(sum(output[[1]])/sum(output[[2]]), 1)
    expect_equivalent(u1, as.matrix(output[[1]]))

    # Works with blocking factors.
    set.seed(0)
    output <- downsampleBatches(u1, u1*10, u1 * 5,  u1*2, block=c(1,2,1,2))
    expect_equal(colSums(output[[1]])/colSums(output[[3]]), rep(1, ncol(output[[1]])))
    expect_equivalent(u1, as.matrix(output[[1]]))
    expect_equal(colSums(output[[2]])/colSums(output[[4]]), rep(1, ncol(output[[1]])))
    expect_equivalent(u1*2, as.matrix(output[[4]]))

    set.seed(0)
    ref1 <- downsampleBatches(u1, u1*5)
    ref2 <- downsampleBatches(u1*10, u1*2)
    expect_identical(output[c(1,3)], ref1)
    expect_identical(output[c(2,4)], ref2)

    # Checking that it's a no-op when the coverage is the same
    # (aside from the type conversion).
    mat <- as(u1, "dgCMatrix")
    expect_equal(downsampleBatches(u1, u1), List(mat, mat))
})

set.seed(5001)
test_that("downsampling batches gives consistent results with a single object", {
    u1 <- matrix(rpois(20000, 5), ncol=100)
    u2 <- matrix(rpois(40000, 1), ncol=200)
    combined <- cbind(u1, u2)
    batch <- rep(1:2, c(ncol(u1), ncol(u2)))

    for (method in c("mean", "median", "geomean")) {
        set.seed(100)
        output <- downsampleBatches(u1, u2, method=method)
        set.seed(100)
        output2 <- downsampleBatches(combined, batch=batch, method=method)
        expect_identical(output2, do.call(cbind, as.list(output)))
    }

    # Works correctly with bycol=FALSE.
    set.seed(100)
    output <- downsampleBatches(u1, u2, bycol=FALSE)
    set.seed(100)
    output2 <- downsampleBatches(combined, batch=batch, bycol=FALSE)

    # Responds correctly with scrambled ordering.
    scrambler <- sample(batch)
    new.indices <- integer(ncol(combined))
    new.indices[scrambler==1] <- seq_len(ncol(u1))
    new.indices[scrambler==2] <- ncol(u1) + seq_len(ncol(u2))

    set.seed(100)
    ref <- downsampleBatches(combined, batch=batch)
    set.seed(100)
    alt <- downsampleBatches(combined[,new.indices], batch=batch[new.indices])
    expect_identical(ref[,new.indices], alt)

    # Handles the blocking correctly.
    v1 <- matrix(rpois(20000, 10), ncol=100)
    v2 <- matrix(rpois(40000, 2), ncol=200)
    set.seed(100)
    ref <- downsampleBatches(u1, u2, v1, v2, block=c(1,2,1,2))
    set.seed(100)
    mult <- c(ncol(u1), ncol(u2), ncol(v1), ncol(v2))
    alt <- downsampleBatches(cbind(u1, u2, v1, v2), batch=rep(1:4, mult), block=rep(c(1,2,1,2), mult))
    expect_identical(do.call(cbind, as.list(ref)), alt)

    expect_error(downsampleBatches(combined), "must be specified")
    expect_error(downsampleBatches(combined, batch="A"), "must be equal")
    expect_error(downsampleBatches(combined, batch=batch, block="A"), "must be equal")
})

set.seed(5001)
test_that("downsampling batches handles a diverse set of other inputs", {
    u1 <- matrix(rpois(20000, 5), ncol=100)
    u2 <- matrix(rpois(40000, 1), ncol=200)

    set.seed(0)
    ref <- downsampleBatches(u1, u2)

    set.seed(0)
    alt <- downsampleBatches(list(u1, u2))
    expect_identical(ref, alt)

    set.seed(0)
    alt <- downsampleBatches(SummarizedExperiment(u1), SummarizedExperiment(u2))
    expect_identical(ref, alt)
})

