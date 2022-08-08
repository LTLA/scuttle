# tests for cell-based pre-processing functions.
# library(scuttle); library(testthat); source("setup.R"); source("test-sum-across-cells.R")

library(Matrix)
library(DelayedArray)

##########################################################

set.seed(10003)
test_that("we can summarise counts at cell cluster level", {
    ids <- sample(ncol(sce)/2, ncol(sce), replace=TRUE)
    out <- sumCountsAcrossCells(sce, ids)

    expect_identical(assay(out), colsum(counts(sce), ids))
    expect_identical(colnames(out), as.character(sort(unique(ids)))) # numeric ordering is preserved.
    expect_identical(sort(unique(ids)), out$ids)

    out2 <- sumCountsAcrossCells(counts(sce), ids)
    expect_identical(out, out2)

    # Robust to column names.
    copy <- sce
    colnames(copy) <- paste0("CELL", seq_len(ncol(copy)))
    expect_identical(sumCountsAcrossCells(copy, ids), out)

    # assay.type= works correctly.
    alt <- sce
    assayNames(alt) <- "whee"
    out2 <- sumCountsAcrossCells(alt, ids, assay.type="whee")
    expect_identical(out, out2)

    # Respects levels properly.
    fids <- factor(ids, levels=rev(sort(unique(ids))))
    fout <- sumCountsAcrossCells(sce, fids)
    fout <- fout[,ncol(fout):1]
    fout$ids <- as.integer(levels(fout$ids))[fout$ids]
    expect_identical(out, fout)

    # Drops unused levels.
    fids <- factor(ids, levels=1:100)
    fout <- sumCountsAcrossCells(sce, fids)
    expect_identical(fout$ids, sort(unique(fids)))
    fout$ids <- as.integer(fout$ids)
    expect_identical(out, fout)

    # Handles NA's correctly.
    ids2 <- sample(LETTERS, ncol(sce), replace=TRUE)
    out2 <- sumCountsAcrossCells(sce, ids2)

    ids3 <- ids2
    ids3[ids3=="A"] <- NA
    out3 <- sumCountsAcrossCells(sce, ids3)

    expect_identical(out2[,setdiff(colnames(out2), "A")], out3)

    all.na <- ids
    all.na[] <- NA
    out3 <- sumCountsAcrossCells(sce, all.na)
    expect_identical(ncol(out3), 0L)
    expect_identical(rownames(out3), rownames(sce))

    # Doesn't store numbers fi we don't ask.
    out <- sumCountsAcrossCells(sce, ids, store.number=NULL)
    expect_null(out$ncells)
})

set.seed(10004)
test_that("by-cell count summarization behaves with other classes", {
    ids <- sample(ncol(sce)/2, ncol(sce), replace=TRUE)
    ref <- sumCountsAcrossCells(sce, ids)

    # Handles sparse matrices properly.
    sparsified <- sce
    counts(sparsified) <- as(counts(sparsified), "dgCMatrix")
    spack <- sumCountsAcrossCells(sparsified, ids)
    expect_identical(ref, spack)

    unknown <- sce
    counts(unknown) <- as(counts(unknown), "dgTMatrix")
    spack <- sumCountsAcrossCells(unknown, ids)
    expect_identical(ref, spack)

    # Handles DelayedArrays properly.
    delayed <- sce
    counts(delayed) <- DelayedArray(counts(delayed))
    dack <- sumCountsAcrossCells(delayed, ids)
    expect_equivalent(ref, dack)

    # Handles _sparse_ DelayedArrays properly.
    sdelayed <- sce
    counts(sdelayed) <- DelayedArray(as(counts(sce), "dgCMatrix"))
    sdack <- sumCountsAcrossCells(sdelayed, ids)
    expect_equivalent(ref, sdack)
})

set.seed(100041)
test_that("by-cell count summarization handles parallelization properly", {
    ids <- sample(ncol(sce)/2, ncol(sce), replace=TRUE)
    ref <- sumCountsAcrossCells(sce, ids)
    
    alt <- sumCountsAcrossCells(sce, ids, BPPARAM=safeBPParam(2))
    expect_identical(alt, ref)

    alt <- sumCountsAcrossCells(sce, ids, BPPARAM=safeBPParam(3))
    expect_identical(alt, ref)
})

set.seed(10004001)
test_that("by-cell count summarization behaves with subsetting", {
    ids <- sample(LETTERS[1:5], ncol(sce), replace=TRUE)

    expect_identical(sumCountsAcrossCells(counts(sce), ids, subset.row=10:1),
        sumCountsAcrossCells(counts(sce), ids)[10:1,])

    expect_identical(sumCountsAcrossCells(counts(sce), ids, subset.col=2:15),
        sumCountsAcrossCells(counts(sce)[,2:15], ids[2:15]))
})

set.seed(100040002)
test_that("by-cell count summarization works with various average types", {
    ids <- sample(LETTERS[1:5], ncol(sce), replace=TRUE)
    ref <- sumCountsAcrossCells(sce, ids)

    # Handles vanilla averaging:
    out2 <- sumCountsAcrossCells(sce, ids, average=TRUE)
    expect_identical(assay(out2), t(t(colsum(counts(sce), ids))/as.integer(table(ids))))
    expect_identical(colData(out2), colData(ref))

    out3 <- sumCountsAcrossCells(sce, ids, average="mean")
    expect_identical(out2, out3)

    # Handles medianizing.
    out4 <- sumCountsAcrossCells(sce, ids, average="median")
    expect_false(identical(out2, out4))

    # 'none' is the same as FALSE.
    expect_identical(ref, sumCountsAcrossCells(sce, ids, average="none"))
})

set.seed(1000401)
test_that("Aggregation across cells works correctly with DFs", {
    # One factor.
    ids <- sample(ncol(sce)/2, ncol(sce), replace=TRUE)
    ref <- sumCountsAcrossCells(sce, ids)
    out <- sumCountsAcrossCells(sce, DataFrame(X=ids))

    expect_identical(colnames(ref), as.character(out$X))
    expect_equivalent(assay(ref), assay(out))
    expect_identical(out$ncells, ref$ncells)
    expect_identical(out$ncells, as.integer(table(ids)))

    # Two factors.
    extra <- sample(LETTERS[1:3], ncol(sce), replace=TRUE)
    combined <- paste0(ids, "-", extra)
    ref <- sumCountsAcrossCells(sce, combined)
    df <- DataFrame(X=ids, Y=extra)
    out <- sumCountsAcrossCells(sce, df)

    post.combined <- paste0(out$X, "-", out$Y)
    expect_identical(sort(colnames(ref)), sort(post.combined))
    m <- match(colnames(ref), post.combined)
    expect_equivalent(assay(ref), assay(out)[,m])

    expect_identical(order(colData(out)), seq_len(ncol(out))) # output is ordered.
    expect_identical(out$ncells, as.integer(table(selfmatch(sort(df)))))

    ref <- sumCountsAcrossCells(sce, combined, average=TRUE)
    out <- sumCountsAcrossCells(sce, df, average=TRUE)
    expect_equivalent(assay(ref), assay(out)[,m])

    # Handles NAs correctly.
    extra[1] <- NA
    ids[2] <- NA
    df <- DataFrame(X=ids, Y=extra)

    ref <- sumCountsAcrossCells(sce[,-(1:2)], df[-(1:2),])
    out <- sumCountsAcrossCells(sce, df)
    expect_equal(assay(ref), assay(out))
    expect_equal(colData(ref), colData(out))

    out2 <- sumCountsAcrossCells(sce, df, subset.col=-(1:2))
    expect_equal(out, out2)
})

##########################################################

set.seed(100041)
test_that("Aggregation across cells works correctly for SCEs", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    alt <- aggregateAcrossCells(sce, ids)
    expect_identical(colnames(alt), sort(unique(ids)))

    ref <- sumCountsAcrossCells(counts(sce), ids)
    expect_identical(counts(alt), assay(ref))
    expect_identical(alt$ncells, ref$ncells)
    expect_identical(alt$ids, ref$ids)

    # Behaves in the presence of multiple assays.
    normcounts(sce) <- normalizeCounts(sce, log=FALSE)
    alt2 <- aggregateAcrossCells(sce, ids)
    expect_identical(alt, alt2)

    sce <- logNormCounts(sce, log=FALSE)
    alt3 <- aggregateAcrossCells(sce, ids, use.assay.type=c("counts", "normcounts"))
    expect_identical(counts(alt), counts(alt3))

    ref <- sumCountsAcrossCells(sce, ids, assay.type="normcounts")
    expect_identical(normcounts(alt3), assay(ref))

    # Works when the count matrix is not the first.
    sce2 <- sce
    assays(sce2) <- assays(sce2)[2:1]
    alt4 <- aggregateAcrossCells(sce2, ids)
    ref <- sumCountsAcrossCells(counts(sce), ids)
    expect_identical(counts(alt4), assay(ref))
})

set.seed(1000401)
test_that("Aggregation across cells works correctly with altExps", {
    ids <- paste0("CLUSTER_", sample(10, ncol(sce), replace=TRUE))
    copy <- sce
    altExp(copy, "THING") <- sce
    counts(altExp(copy)) <- counts(altExp(copy)) * 2

    # Doesn't pass along by default.
    suppressWarnings(agg <- aggregateAcrossCells(copy, ids))
    expect_identical(altExpNames(agg), character(0))

    # But apply'ing works correctly.
    agg0 <- applySCE(copy, aggregateAcrossCells, ids=ids)
    expect_identical(counts(agg0), counts(agg))
    expect_identical(counts(altExp(agg0, "THING")), counts(agg)*2)

    # Other options work correctly.
    agg4 <- applySCE(copy, FUN=aggregateAcrossCells, ids=ids, WHICH=1, use.altexps=NULL)
    expect_identical(altExpNames(agg4), "THING")
    agg5 <- applySCE(copy, FUN=aggregateAcrossCells, ids=ids, WHICH="THING", use.altexps=NULL)
    expect_identical(altExpNames(agg5), "THING")
})

set.seed(1000401)
test_that("Aggregation across cells works correctly with reducedDims", {
    ids <- paste0("CLUSTER_", sample(20, ncol(sce), replace=TRUE))
    copy <- sce
    reducedDim(copy, "PCA") <- t(assay(sce)[1:3,])
    reducedDim(copy, "TSNE") <- t(assay(sce)[1:10,])

    # Responds to the average settings.
    agg <- aggregateAcrossCells(copy, ids, statistics="mean")
    expect_identical(reducedDim(agg, "PCA"), t(assay(agg)[1:3,]))
    expect_identical(reducedDim(agg, "TSNE"), t(assay(agg)[1:10,]))

    agg2 <- aggregateAcrossCells(copy, ids, statistics="sum")
    expect_identical(reducedDims(agg2), reducedDims(agg))

    agg3 <- aggregateAcrossCells(copy, ids, statistics="median", dimred.stats="median")
    expect_identical(reducedDim(agg3, "PCA"), t(assay(agg3)[1:3,]))
    expect_identical(reducedDim(agg3, "TSNE"), t(assay(agg3)[1:10,]))
    expect_false(identical(agg3, agg))

    # Behaves with NAs.
    ids2 <- ids
    failed <- ids2==ids2[1]
    ids2[failed] <- NA
    expect_identical(
        aggregateAcrossCells(copy, ids2, statistics="mean"),
        aggregateAcrossCells(copy[,!failed], ids[!failed], statistics="mean")
    )

    # Other options work correctly.
    agg1 <- aggregateAcrossCells(copy, ids, use.dimred=1)
    expect_identical(reducedDimNames(agg1), "PCA")
    expect_error(aggregateAcrossCells(copy, ids, use.dimred=10), 'use.dimred')

    agg2 <- aggregateAcrossCells(copy, ids, use.dimred="TSNE")
    expect_identical(reducedDimNames(agg2), "TSNE")
    expect_error(aggregateAcrossCells(copy, ids, use.dimred="WHEE"), 'use.dimred')

    agg0 <- aggregateAcrossCells(sce, ids, use.dimred=FALSE)
    expect_identical(reducedDimNames(agg0), character(0))
})

set.seed(1000411)
test_that("Aggregation across cells works correctly for SCEs with DFs", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    extra <- sample(LETTERS[1:3], ncol(sce), replace=TRUE)

    combined <- DataFrame(X=ids, Y=extra)
    agg <- aggregateAcrossCells(sce, combined)
    ref <- sumCountsAcrossCells(counts(sce), combined)

    expect_identical(counts(agg), assay(ref))
    expect_identical(agg$X, ref$X)
    expect_identical(agg$Y, ref$Y)
    expect_identical(agg$ncells, ref$ncells)

    # Same for alternative experiments.
    copy <- sce
    altExp(copy, "THING") <- sce
    counts(altExp(copy)) <- counts(altExp(copy)) * 2

    agg <- applySCE(copy, aggregateAcrossCells, ids=combined)
    expect_identical(counts(agg), assay(ref))
    expect_identical(counts(altExp(agg, "THING")), assay(ref)*2)
    expect_identical(ref$X, altExp(agg)$X)
    expect_identical(ref$Y, altExp(agg)$Y)

    expect_identical(agg$ncells, ref$ncells)
    expect_identical(altExp(agg)$ncells, ref$ncells)

    # Same for reduced dimensions.
    copy <- sce
    reducedDim(copy, "PCA") <- t(assay(sce)[1:3,])
    agg <- aggregateAcrossCells(copy, combined, statistics="mean")
    expect_identical(reducedDim(agg, "PCA"), t(assay(agg)[1:3,]))
})

set.seed(1000412)
test_that("Aggregation across cells works correctly with custom coldata acquisition", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    sce$thing <- seq_len(ncol(sce))

    # Defaults to partial NA's.
    alt <- aggregateAcrossCells(sce, ids)
    tab <- table(ids)
    expect_identical(colnames(alt)[is.na(alt$thing)], names(tab)[tab > 1])
    expect_equivalent(alt$thing[!is.na(alt$thing)], sce$thing[match(names(tab)[tab==1], ids)])

    # Defaults to a sensible value if we enforce identity within each group.
    sce$thing2 <- ids
    alt <- aggregateAcrossCells(sce, ids)
    expect_equivalent(colnames(alt), alt$thing2)

    # Responds to taking the first.
    alt <- aggregateAcrossCells(sce, ids, coldata_merge=function(x) head(x, 1))
    expect_equivalent(alt$thing, as.integer(by(sce$thing, ids, head, n=1)))
    expect_equivalent(alt$Mutation_Status, as.character(
        by(data.frame(sce$Mutation_Status, stringsAsFactors=FALSE), ids, FUN=head, n=1))
    )
    expect_identical(colnames(alt), sort(unique(ids)))

    # Responds to taking the sum.
    alt <- aggregateAcrossCells(sce, ids, coldata_merge=list(thing=sum))
    expect_equivalent(alt$thing, as.integer(by(sce$thing, ids, sum)))
    expect_identical(colnames(alt), sort(unique(ids)))

    alt <- aggregateAcrossCells(sce, ids, coldata_merge=list(Cell_Cycle=function(x) paste(x, collapse="")))
    expect_type(alt$Cell_Cycle, "character")

    # Setting FALSE works corectly.
    alt <- aggregateAcrossCells(sce, ids, coldata_merge=FALSE)
    expect_identical(colnames(colData(alt)), c("ids", "ncells"))
    alt <- aggregateAcrossCells(sce, ids, coldata_merge=list(thing=FALSE))
    expect_identical(alt$thing, NULL)
})

set.seed(1000412)
test_that("Aggregation across cells works correctly with different colData types", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))

    # Handles factors correctly.
    alt <- sce
    alt$Cell_Cycle <- factor(alt$Cell_Cycle)
    alt <- aggregateAcrossCells(alt, factor(ids))
    expect_true(is.factor(alt$Cell_Cycle))
    expect_true(is.factor(alt$ids))

    alt <- sce
    alt$Cell_Cycle <- factor(alt$Cell_Cycle)
    alt <- aggregateAcrossCells(alt, rep(NA, ncol(alt)))
    expect_true(is.factor(alt$Cell_Cycle))

    # Handles nested DFs correctly.
    alt <- sce
    alt$nested <- DataFrame(ID=ids)
    alt <- aggregateAcrossCells(alt, ids)
    expect_s4_class(alt$nested, "DFrame")
    expect_identical(alt$nested$ID, alt$ids)

    alt <- sce
    alt$nested <- DataFrame(ID=runif(ncol(alt)))
    alt <- aggregateAcrossCells(alt, ids)
    expect_s4_class(alt$nested, "DFrame")
    expect_type(alt$nested$ID, "double")
    expect_true(any(is.na(alt$nested$ID)))

    alt <- sce
    alt$nested <- DataFrame(ID=runif(ncol(alt)))
    alt <- aggregateAcrossCells(alt, rep(NA, ncol(alt)))
    expect_s4_class(alt$nested, "DFrame")
    expect_type(alt$nested$ID, "double")
})

set.seed(10004121)
test_that("Aggregation across cells works correctly with suffixing", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    alt <- aggregateAcrossCells(sce, ids, suffix=TRUE)
    expect_identical(assayNames(alt), "counts.sum")

    sce <- logNormCounts(sce)
    alt <- aggregateAcrossCells(sce, ids, statistics=c("mean", "median"), use.assay.type=TRUE)
    expect_identical(assayNames(alt), c("counts.mean", "counts.median", "logcounts.mean", "logcounts.median"))

    # Works for reddims.
    reducedDims(sce) <- list(PCA=matrix(runif(ncol(sce)*2), ncol=2), TSNE=matrix(rnorm(ncol(sce)*2), ncol=2))
    alt <- aggregateAcrossCells(sce, ids)
    expect_identical(reducedDimNames(alt), reducedDimNames(sce))

    alt <- aggregateAcrossCells(sce, ids, suffix=TRUE)
    expect_identical(reducedDimNames(alt), paste0(reducedDimNames(sce), ".mean"))
    expect_identical(assayNames(alt), "counts.sum") # passed along.

    alt <- aggregateAcrossCells(sce, ids, use.dimred="PCA", dimred.stats=c("mean", "median"))
    expect_identical(reducedDimNames(alt), c("PCA.mean", "PCA.median"))
})

set.seed(1000413)
test_that("Aggregation across cells doesn't choke on syntactically invalid names", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    sce$`batman and robin` <- 2
    alt <- aggregateAcrossCells(sce, ids)
    expect_identical(sum(alt$`batman and robin`==2), ncol(alt))
})

set.seed(100042)
test_that("Aggregation across cells works correctly for SEs", {
    ids <- paste0("CLUSTER_", sample(ncol(sce)/2, ncol(sce), replace=TRUE))
    alt <- aggregateAcrossCells(sce, ids)
    expect_identical(colnames(alt), sort(unique(ids)))

    ref <- sumCountsAcrossCells(counts(sce), ids)
    expect_identical(counts(alt), assay(ref))
    expect_identical(alt$ids, ref$ids)
    expect_identical(alt$ncells, ref$ncells)
})

##########################################################

test_that("numDetectedAcrossCells works as expected", {
    ids <- sample(LETTERS[1:5], ncol(sce), replace=TRUE)

    out <- numDetectedAcrossCells(counts(sce), ids)
    expect_equal(assay(out), colsum((counts(sce) > 0)+0, ids))
    expect_identical(out$ids, colnames(out))
    out <- numDetectedAcrossCells(counts(sce), ids, average=TRUE)
    expect_identical(assay(out), t(t(colsum((counts(sce) > 0)+0, ids))/as.integer(table(ids))))

    # Checking that it works direclty with SCEs.
    expect_equal(numDetectedAcrossCells(counts(sce), ids),
        numDetectedAcrossCells(sce, ids))
    expect_equal(numDetectedAcrossCells(counts(sce), ids, average=TRUE),
        numDetectedAcrossCells(sce, ids, average=TRUE))

    # Checking that subsetting works.
    expect_identical(numDetectedAcrossCells(counts(sce), ids, subset.row=10:1),
        numDetectedAcrossCells(counts(sce), ids)[10:1,])

    expect_identical(numDetectedAcrossCells(counts(sce), ids, subset.col=2:15),
        numDetectedAcrossCells(counts(sce)[,2:15], ids[2:15]))

    ids[c(1,3,5,6)] <- NA
    expect_identical(numDetectedAcrossCells(counts(sce), ids),
        numDetectedAcrossCells(counts(sce)[,!is.na(ids)], ids[!is.na(ids)]))

    # Comparing to sumCountsAcrossCells.
    expect_equal(numDetectedAcrossCells(counts(sce), ids),
        sumCountsAcrossCells((counts(sce) > 0)+0, ids))
    expect_equal(numDetectedAcrossCells(counts(sce), ids, average=TRUE),
        sumCountsAcrossCells((counts(sce) > 0)+0, ids, average=TRUE))
})

test_that("numDetectedAcrossCells handles other matrix classes", {
    thing <- matrix(rpois(2000, lambda=0.5), ncol=100, nrow=20)
    ids <- sample(LETTERS[1:6], ncol(thing), replace=TRUE)

    ref <- numDetectedAcrossCells(thing, ids)
    expect_equal(rowSums(assay(ref)), rowSums(thing > 0)) # basic sanity check.

    sparse <- as(thing, 'dgCMatrix')
    expect_equal(numDetectedAcrossCells(sparse, ids), ref)

    delayed <- DelayedArray(thing)
    expect_equal(numDetectedAcrossCells(delayed, ids), ref)
})
