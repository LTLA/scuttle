# This tests the whichNonZero utility function.
# library(scuttle); library(testthat); source("test-which-nonzero.R")

library(DelayedArray)

test_that("whichNonZero works as expected for all matrix types", {
    stuff <- Matrix::rsparsematrix(1000, 1000, density=0.01)
    ref <- whichNonZero(as.matrix(stuff))

    out <- whichNonZero(stuff)
    expect_identical(ref, out)
    
    out <- whichNonZero(as(stuff, "dgTMatrix"))
    expect_identical(ref, out)

    out <- whichNonZero(DelayedArray(stuff))
    expect_identical(ref, out)
})

expect_identical_sorted <- function(left, right) {
    o <- order(left$i, left$j)
    left <- lapply(left, "[", i=o)
    o <- order(right$i, right$j)
    right <- lapply(right, "[", i=o)
    expect_identical(left, right)
}

test_that("whichNonZero behaves correctly for DelayedArrays", {
    basic <- matrix(rpois(1000, 0.1), ncol=10)
    ref <- whichNonZero(basic)

    # Again, a basic check.
    thing <- DelayedArray(basic)
    out <- whichNonZero(thing)
    expect_identical_sorted(ref, out)

    # Constraining the grid.
    setAutoBlockSize(ncol(thing) * 4L)
    out <- whichNonZero(thing)
    expect_identical_sorted(ref, out)

    # Checking the other dimension of the viewports.
    alt <- DelayedArray(t(basic))
    setAutoBlockSize(nrow(alt) * 4L)
    out <- whichNonZero(alt)
    expect_identical_sorted(out, whichNonZero(t(basic)))

    # Parallelizing.
    out <- whichNonZero(thing, BPPARAM=BiocParallel::SnowParam(3))
    expect_identical_sorted(ref, out)

    setAutoBlockSize() # resetting.
})
