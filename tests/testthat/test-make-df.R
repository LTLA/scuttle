## Tests for the visualization variable picker.
## library(scuttle); library(testthat); source("setup.R"); source("test-make-df.R")

example_sce <- normed

altExp(example_sce, "thing") <- normed[1:10,]
rownames(altExp(example_sce)) <- paste0(rownames(altExp(example_sce)), ".0")

altExp(example_sce, "other") <- normed[1:5,]
rownames(altExp(example_sce, 2)) <- paste0(rownames(altExp(example_sce, 2)), "-R")

set.seed(1312313)
rowData(example_sce) <- DataFrame(HAPPY=runif(nrow(example_sce)), SAD=rbinom(nrow(example_sce), 1, 0.5)==1)

test_that("makePerCellDF works as expected", {
    df1 <- makePerCellDF(example_sce, features="Gene_0001")
    expect_identical(df1$Mutation_Status, example_sce$Mutation_Status)
    expect_identical(df1$Gene_0001, unname(logcounts(example_sce)["Gene_0001",]))

    # Works with reduced dimensions and size factors.
    reducedDim(example_sce, "PCA") <- matrix(runif(ncol(example_sce) * 5), ncol=5)
    df3 <- makePerCellDF(example_sce)
    expect_identical(df3$sizeFactor, unname( sizeFactors(example_sce)))
    expect_identical(df3$PCA.1, unname(reducedDim(example_sce)[,1]))
    expect_identical(df3$PCA.2, unname(reducedDim(example_sce)[,2]))

    df3b <- makePerCellDF(example_sce, use_dimred=FALSE)
    expect_true(all(!grepl("PCA", colnames(df3b))))

    # Works with alternative experiments.
    df4 <- makePerCellDF(example_sce, features="Gene_0001.0", use.altexps=TRUE, prefix.altexps=TRUE)
    expect_identical(df4$thing.Gene_0001, unname(logcounts(altExp(example_sce))["Gene_0001.0",]))

    df4b <- makePerCellDF(example_sce, use.altexps="other", prefix.altexps=TRUE)
    expect_true(all(!grepl("^thing\\.", colnames(df4b))))

    # Handles edge cases gracefully.
    stripped <- example_sce
    colData(stripped) <- NULL
    reducedDims(stripped) <- NULL

    df0 <- makePerCellDF(stripped)
    expect_identical(ncol(df0), 0L)
    expect_identical(rownames(df0), colnames(example_sce))

    # Do not fix my names unless requested.
    rownames(example_sce) <- paste0("+", seq_len(nrow(example_sce)))
    df0b <- makePerCellDF(example_sce, features=rownames(example_sce)[1:10])
    expect_identical(tail(colnames(df0b), 10), rownames(example_sce)[1:10])

    reducedDimNames(example_sce)[1] <- "+-PCA"
    df0c <- makePerCellDF(example_sce)
    expect_true("+-PCA.1" %in% colnames(df0c))

    df0d <- makePerCellDF(example_sce, features=rownames(example_sce)[1:10], check.names=TRUE)
    expect_true("X..PCA.1" %in% colnames(df0d))
    expect_identical(tail(colnames(df0d), 10), make.names(rownames(example_sce)[1:10]))
})

test_that("makePerCellDF works with rowname swapping", {
    rowData(example_sce)$alias <- sprintf("FEATURE_%s", seq_len(nrow(example_sce)))
    df1 <- makePerCellDF(example_sce, features="FEATURE_1", swap.rownames="alias")
    expect_identical(df1$FEATURE_1, unname(logcounts(example_sce)[1,]))

    # Combines easily with altexps.
    rowData(altExp(example_sce))$alias <- sprintf("thing_%s", seq_len(nrow(altExp(example_sce))))
    df1 <- makePerCellDF(example_sce, features=c("FEATURE_10", "thing_9"), swap.rownames="alias")
    expect_identical(df1$FEATURE_10, unname(logcounts(example_sce)[10,]))
    expect_identical(df1$thing_9, unname(logcounts(altExp(example_sce))[9,]))

    df1 <- makePerCellDF(example_sce, features=c("FEATURE_10", "thing_9"), swap.rownames="alias", prefix.altexps=TRUE)
    expect_identical(df1$FEATURE_10, unname(logcounts(example_sce)[10,]))
    expect_identical(df1$thing.thing_9, unname(logcounts(altExp(example_sce))[9,]))

    # Graceully handles absent swaps.
    df1 <- makePerCellDF(example_sce, features="FEATURE_1", swap.rownames="alchemist")
    expect_false(any(grepl("gene|feature", colnames(df1), ignore.case=TRUE)))
})

test_that("makePerFeatureDF works as expected", {
    rowData(example_sce)$foo <- runif(nrow(example_sce))
    rowData(example_sce)$bar <- sample(LETTERS, nrow(example_sce), replace=TRUE)

    df1 <- makePerFeatureDF(example_sce, cells="Cell_001")
    expect_identical(df1$foo, rowData(example_sce)$foo)
    expect_identical(df1$bar, rowData(example_sce)$bar)
    expect_identical(df1$Cell_001, unname(logcounts(example_sce)[,"Cell_001"]))

    # Handles edge cases gracefully.
    stripped <- example_sce
    rowData(stripped) <- NULL
    df0 <- makePerFeatureDF(stripped)
    expect_identical(ncol(df0), 0L)
    expect_identical(rownames(df0), rownames(example_sce))

    # Do not fix my names.
    colnames(example_sce) <- paste0("+", seq_len(ncol(example_sce)))
    df0b <- makePerFeatureDF(example_sce, cells=colnames(example_sce)[1:10])
    expect_identical(colnames(df0b)[1:10], colnames(example_sce)[1:10])

    df0c <- makePerFeatureDF(example_sce, cells=colnames(example_sce)[1:10], check.names=TRUE)
    expect_identical(colnames(df0c)[1:10], make.names(colnames(example_sce)[1:10]))
})

test_that("makePer*DF functions work for non-ordinary matrices", {
    logcounts(example_sce) <- as(logcounts(example_sce), "dgCMatrix")

    df1 <- makePerCellDF(example_sce, features=c("Gene_0001", "Gene_0010", "Gene_0100"))
    expect_identical(df1$Gene_0001, unname(logcounts(example_sce)["Gene_0001",]))
    expect_identical(df1$Gene_0010, unname(logcounts(example_sce)["Gene_0010",]))
    expect_identical(df1$Gene_0100, unname(logcounts(example_sce)["Gene_0100",]))

    df1 <- makePerFeatureDF(example_sce, cells=c("Cell_001", "Cell_010", "Cell_100"))
    expect_identical(df1$Cell_001, unname(logcounts(example_sce)[,"Cell_001",]))
    expect_identical(df1$Cell_010, unname(logcounts(example_sce)[,"Cell_010",]))
    expect_identical(df1$Cell_100, unname(logcounts(example_sce)[,"Cell_100",]))
})
