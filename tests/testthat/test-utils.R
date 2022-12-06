# This tests the various internal utilities.
# library(testthat); library(scuttle); source("setup.R"); source("test-utils.R")

test_that("converting subset vectors to indices is correct", {
    M <- matrix(rnorm(1000), ncol=20)
    colnames(M) <- head(LETTERS, ncol(M))
    rownames(M) <- sprintf("GENE_%s", seq_len(nrow(M)))

    # For matrix rows:
    expect_identical(.subset2index(NULL, M, byrow=TRUE), seq_len(nrow(M)))
    expect_identical(.subset2index(1:5, M, byrow=TRUE), 1:5)
    expect_identical(.subset2index(rownames(M)[10:5], M, byrow=TRUE), 10:5)
    expect_identical(.subset2index(seq_len(nrow(M)) %in% 20:50, M, byrow=TRUE), 20:50)
    expect_error(.subset2index(1000:2000, M, byrow=TRUE), "invalid")

    # For matrix columns:
    expect_identical(.subset2index(NULL, M, byrow=FALSE), seq_len(ncol(M)))
    expect_identical(.subset2index(1:5, M, byrow=FALSE), 1:5)
    expect_identical(.subset2index(colnames(M)[10:5], M, byrow=FALSE), 10:5)
    expect_identical(.subset2index(seq_len(ncol(M)) %in% 10:15, M, byrow=FALSE), 10:15)
    expect_error(.subset2index(1000:2000, M, byrow=FALSE), "invalid")

    # For vectors:
    V <- M[,1]
    expect_identical(.subset2index(NULL, V, byrow=NA), seq_along(V))
    expect_identical(.subset2index(1:5, V, byrow=NA), 1:5)
    expect_identical(.subset2index(names(V)[10:5], V, byrow=NA), 10:5)
    expect_identical(.subset2index(seq_along(V) %in% 10:15, V, byrow=NA), 10:15)
    expect_error(.subset2index(1000:2000, V, byrow=NA), "invalid")

    # Factors:
    expect_identical(.subset2index(factor(rownames(M)[4:5]), M), 4:5)
    expect_identical(.subset2index(factor(colnames(M)[10:5]), M, byrow=FALSE), 10:5)

    # Edge cases.
    expect_identical(.subset2index(numeric(0), M, byrow=FALSE), integer(0))
    expect_identical(.subset2index(character(0), M, byrow=FALSE), integer(0))
    expect_identical(.subset2index(logical(0), M, byrow=FALSE), integer(0))
})

test_that("job assignment to workers is correct", {
    for (nc in c(3, 6, 11)) {
        P <- safeBPParam(nc)
        n <- bpnworkers(P)

        out <- .assignIndicesToWorkers(100, P)
        expect_identical(length(out), n)
        expect_true(all(lengths(out) >= floor(100/n)))
        expect_identical(unlist(out), seq_len(100))
    }

    # Works with subsetting.
    chosen <- sample(100, 50)
    P <- safeBPParam(11)
    n <- bpnworkers(P)
    out <- .assignIndicesToWorkers(NULL, P, subset=chosen)
    expect_identical(length(out), n)
    expect_true(all(lengths(out) >= floor(length(chosen)/n)))
    expect_identical(unlist(out), chosen)

    chosen <- sample(LETTERS)
    P <- safeBPParam(5)
    n <- bpnworkers(P)
    out <- .assignIndicesToWorkers(NULL, P, subset=chosen)
    expect_identical(length(out), n)
    expect_true(all(lengths(out) >= floor(length(chosen)/n)))
    expect_identical(unlist(out), chosen)

    chosen <- rbinom(25, 1, 0.5)==1
    P <- safeBPParam(3)
    n <- bpnworkers(P)
    out <- .assignIndicesToWorkers(NULL, P, subset=chosen)
    expect_identical(length(out), n)
    expect_true(all(lengths(out) >= floor(sum(chosen)/n)))
    expect_identical(unlist(out), which(chosen))
})


test_that("splitting a vector to workers is correct", {
    X <- runif(99)

    P <- safeBPParam(7)
    out <- .splitVectorByWorkers(X, P)
    expect_identical(unlist(out), X)
    expect_identical(length(out), bpnworkers(P))

    out <- .splitVectorByWorkers(X, safeBPParam(1))
    expect_identical(out[[1]], X)

    # Behaves with subsetting of all flavors.
    i <- 1:50
    expect_identical(
        .splitVectorByWorkers(X, safeBPParam(7), subset=i),
        .splitVectorByWorkers(X[i], safeBPParam(7))
    )

    expect_identical(
        .splitVectorByWorkers(X, safeBPParam(1), subset=i),
        .splitVectorByWorkers(X[i], safeBPParam(1))
    )

    j <- rbinom(99, 1, 0.2)==1
    expect_identical(
        .splitVectorByWorkers(X, safeBPParam(1), subset=j),
        .splitVectorByWorkers(X[j], safeBPParam(1))
    )
})

test_that("splitting a matrix by row is correct", {
    M <- matrix(rnorm(1000), ncol=10)

    P <- safeBPParam(7)
    out <- .splitRowsByWorkers(M, P)
    expect_identical(do.call(rbind, out), M)
    expect_identical(length(out), bpnworkers(P))

    P <- safeBPParam(1)
    out <- .splitRowsByWorkers(M, P)
    expect_identical(out[[1]], M)
    expect_identical(length(out), bpnworkers(P))

    # Behaves with subsetting by row:
    i <- 1:50
    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(7), subset.row=i),
        .splitRowsByWorkers(M[i,], safeBPParam(7))
    )

    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(1), subset.row=i),
        .splitRowsByWorkers(M[i,], safeBPParam(1))
    )

    # Behaves with column subsetting:
    i <- 5:1
    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(7), subset.col=i),
        .splitRowsByWorkers(M[,i], safeBPParam(7))
    )

    expect_identical(
        .splitRowsByWorkers(M, safeBPParam(1), subset.col=i),
        .splitRowsByWorkers(M[,i], safeBPParam(1))
    )
})

test_that("splitting a matrix by column is correct", {
    M <- matrix(rnorm(1000), nrow=20)

    P <- safeBPParam(8)
    out <- .splitColsByWorkers(M, P)
    expect_identical(do.call(cbind, out), M)
    expect_identical(length(out), bpnworkers(P))

    P <- safeBPParam(1)
    out <- .splitColsByWorkers(M, P)
    expect_identical(out[[1]], M)
    expect_identical(length(out), bpnworkers(P))

    # Behaves with subsetting by row:
    i <- 1:10
    expect_identical(
        .splitColsByWorkers(M, safeBPParam(4), subset.row=i),
        .splitColsByWorkers(M[i,], safeBPParam(4))
    )

    expect_identical(
        .splitColsByWorkers(M, safeBPParam(1), subset.row=i),
        .splitColsByWorkers(M[i,], safeBPParam(1))
    )

    # Behaves with column subsetting:
    i <- sample(ncol(M), 20)
    expect_identical(
        .splitColsByWorkers(M, safeBPParam(4), subset.col=i),
        .splitColsByWorkers(M[,i], safeBPParam(4))
    )

    expect_identical(
        .splitColsByWorkers(M, safeBPParam(1), subset.col=i),
        .splitColsByWorkers(M[,i], safeBPParam(1))
    )
})

test_that("unpacking lists is correct", {
    expect_identical(.unpackLists(1, 2, 3), list(1,2,3))

    expect_identical(.unpackLists(1, list(2, 3)), list(1,2,3))

    expect_identical(.unpackLists(list(1, 2, 3)), list(1,2,3))

    # Works with the DataFrame special case.
    df1 <- DataFrame(X=1)
    df2 <- DataFrame(Y=2)
    df3 <- DataFrame(Z=2)

    expect_identical(.unpackLists(df1, df2, df3), list(df1,df2,df3))

    expect_identical(.unpackLists(df1, list(df2, df3)), list(df1,df2,df3))

    expect_identical(.unpackLists(list(df1, df2, df3)), list(df1,df2,df3))

    # Works with names.
    expect_identical(.unpackLists(A=1, B=2, C=3), list(A=1, B=2, C=3))

    expect_identical(.unpackLists(A=1, list(B=2, C=3)), list(A=1, B=2, C=3))

    expect_identical(.unpackLists(list(A=1, B=2, C=3)), list(A=1, B=2, C=3))
})
