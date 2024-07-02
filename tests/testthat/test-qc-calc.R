# Test functions for QC calculation.
# library(scuttle); library(testthat); source("setup.R"); source("test-qc-calc.R")

original <- sce

test_that("we can compute standard per-cell QC metrics", {
    df <- perCellQCMetrics(original, flatten=FALSE)
    expect_identical(rownames(df), colnames(original))

    # Testing total metrics for cells.
    expect_equal(df$sum, unname(colSums(counts(original))))
    expect_equal(df$detected, unname(colSums(counts(original) > 0)))

    # Testing percentage metrics for cells.
    N <- c(50, 100, 200, 500)
    df <- perCellQCMetrics(original, percent.top=N, flatten=FALSE)

    for (i in seq_len(ncol(original))) { 
        cur_counts <- counts(original)[,i]
        o <- order(cur_counts, decreasing=TRUE)
        lib_size <- sum(cur_counts)

        for (x in N) { 
            chosen <- o[seq_len(x)]
            expect_equivalent(df$percent.top[i,as.character(x)], sum(cur_counts[chosen])/lib_size * 100) 
        }
    }

    # Flattening works as expected.
    flat <- perCellQCMetrics(original, percent.top=N)
    for (x in N) {
        expect_identical(flat[,paste0("percent.top_", x)], df$percent.top[,as.character(x)])
    }
    expect_identical(rownames(flat), colnames(original))
})

test_that("we can compute standard QC metrics with subsets", {
    ref <- perCellQCMetrics(original, flatten=FALSE)
    df <- perCellQCMetrics(original, subsets = list(set1 = 1:20), flatten=FALSE)
    expect_identical(df[,1:3], ref[,1:3])

    expect_equivalent(df$subsets$set1$sum, colSums(counts(original[1:20,])))
    expect_equivalent(df$subsets$set1$detected, colSums(counts(original[1:20,])> 0))
    expect_equivalent(df$subsets$set1$percent, df$subsets$set1$sum/df$sum * 100)

    # Testing behaviour with multiple feature controls.
    multi_controls <- list(controls1 = 1:20, controls2 = rownames(original)[500:1000])
    df2 <- perCellQCMetrics(original, subsets = multi_controls, flatten=FALSE)

    expect_equivalent(df$subsets$set1, df2$subsets$controls1)
    expect_equivalent(df2$subsets$controls2$sum, colSums(counts(original[500:1000,])))
    expect_equivalent(df2$subsets$controls2$detected, colSums(counts(original[500:1000,])> 0))
    expect_equivalent(df2$subsets$controls2$percent, df2$subsets$controls2$sum/df2$sum * 100)

    # Flattening works as expected.
    flat <- perCellQCMetrics(original, subsets = list(set1=1:20))
    expect_identical(flat$subsets_set1_sum, df$subsets$set1$sum)
})

test_that("perCellQCMetrics works with alternative experiments", {
    sce <- original
    altExp(sce, "alpha") <- original[1:10,]
    altExp(sce, "bravo") <- original[10:20,]

    ref <- perCellQCMetrics(original, flatten=FALSE)
    df <- perCellQCMetrics(sce, flatten=FALSE)
    expect_identical(df[,1:3], ref[,1:3])

    for (x in altExpNames(sce)) {
        current <- perCellQCMetrics(altExp(sce, x), flatten=FALSE)
        expect_identical(df$altexps[[x]]$sum, current$sum)
        expect_identical(df$altexps[[x]]$detected, current$detected)
        expect_equal(df$altexps[[x]]$percent, current$sum/df$total*100)
    }
    
    expect_identical(df$total, df$sum + df$altexps$alpha$sum + df$altexps$bravo$sum)

    # Ignores experiments that don't have the requested assay.
    copy <- sce
    assayNames(altExp(copy, "alpha")) <- "FOO"
    ignored <- perCellQCMetrics(copy, use.altexps=NULL)
    expect_true(is.null(ignored$altexps_alpha_sum))
    expect_false(is.null(ignored$altexps_bravo_sum))

    # Flattening works as expected.
    flat <- perCellQCMetrics(sce)
    expect_identical(flat$altexps_alpha_sum, df$altexps$alpha$sum)
    expect_identical(flat$altexps_bravo_detected, df$altexps$bravo$detected)

    flat <- perCellQCMetrics(sce, use.altexps="alpha")
    expect_identical(flat$altexps_alpha_sum, df$altexps$alpha$sum)
    expect_null(flat$altexps_beta_sum)
})

test_that("perCellQCMetrics handles silly inputs", {
    expect_error(perCellQCMetrics(original, subsets = list(1:20), flatten=FALSE), "must be named")

    # Doesn't choke with no entries.
    thing <- perCellQCMetrics(original[0,], flatten=FALSE)
    expect_identical(rownames(thing), colnames(original))
    expect_true(all(thing$sum==0L))

    thing2 <- perCellQCMetrics(original[,0], flatten=FALSE)
    expect_identical(nrow(thing2), 0L)
    expect_identical(colnames(thing), colnames(thing2))

    # Percentage holds at the limit.
    df <- perCellQCMetrics(original[1:10,], percent.top=c(50, 100, 200, 500), flatten=FALSE)
    expect_true(all(df$percent.top==100))

    df <- perCellQCMetrics(original, percent.top=integer(0), flatten=FALSE)
    expect_identical(ncol(df$percent.top), 0L)

    # Responds to alternative inputs.
    blah <- sce
    assayNames(blah) <- "whee"
    expect_error(perCellQCMetrics(blah, flatten=FALSE), "counts")
    expect_error(perCellQCMetrics(blah, assay.type="whee", flatten=FALSE), NA)
})

#######################################################################
# Works for per-feature metrics.

test_that("perFeatureQCMetrics works correctly", {
    out <- perFeatureQCMetrics(original, flatten=FALSE)
    expect_identical(rownames(out), rownames(original))
    expect_equal(out$mean, unname(rowMeans(counts(original))))
    expect_equal(out$detected, unname(rowMeans(counts(original) > 0))*100)
})

test_that("we can compute standard QC metrics with cell controls", {
    expect_error(perFeatureQCMetrics(original, subsets = list(1:20), flatten=FALSE), "must be named")

    df <- perFeatureQCMetrics(original, subsets = list(set1 = 1:20), flatten=FALSE)
    sub_counts <- counts(original)[,1:20]

    expect_equal(df$subsets$set1$mean, unname(rowMeans(sub_counts)))
    expect_equal(df$subsets$set1$detected, unname(rowMeans(sub_counts > 0) * 100))
    expect_equal(df$subsets$set1$ratio, df$subsets$set1$mean/df$mean)

    # Testing behaviour with multiple cell controls.
    multi_controls <- list(controls1 = 1:5, controls2 = 10:20)
    df2 <- perFeatureQCMetrics(original, subsets = multi_controls, flatten=FALSE)

    expect_equivalent(df2$subsets$controls2$mean, rowMeans(counts(original[,10:20])))
    expect_equivalent(df2$subsets$controls2$detected, rowMeans(counts(original[,10:20])> 0)*100)
    expect_equivalent(df2$subsets$controls2$ratio, df2$subsets$controls2$mean/df2$mean)

    # Flattening works as expected.
    flat <- perFeatureQCMetrics(original, subsets = list(set1 = 1:20))
    expect_identical(rownames(flat), rownames(original))
    expect_identical(flat$subsets_set1_mean, df$subsets$set1$mean)
    expect_identical(flat$subsets_set1_ratio, df$subsets$set1$ratio)
})

test_that("perFeatureQCmetrics handles silly inputs", {
    expect_error(perFeatureQCMetrics(original, subsets = list(1:20), flatten=FALSE), "must be named")

    # Doesn't choke with no entries.
    thing <- perFeatureQCMetrics(original[,0], flatten=FALSE)
    expect_identical(rownames(thing), rownames(original))
    expect_true(all(thing$sum==0L))

    thing2 <- perFeatureQCMetrics(original[0,], flatten=FALSE)
    expect_identical(nrow(thing2), 0L)
    expect_identical(colnames(thing), colnames(thing2))

    # Responds to alternative inputs.
    blah <- sce
    assayNames(blah) <- "whee"
    expect_error(perFeatureQCMetrics(blah, flatten=FALSE), "counts")
    expect_error(perFeatureQCMetrics(blah, assay.type="whee", flatten=FALSE), NA)
})

#######################################################################
# Responds to special settings: 

test_that("we can compute standard QC metrics on sparse counts matrix", {
    alt <- original

    library(Matrix)
    counts(alt) <- as(counts(alt), "dgCMatrix")

    expect_equal(perCellQCMetrics(alt, flatten=FALSE), perCellQCMetrics(original, flatten=FALSE))
    expect_equal(perFeatureQCMetrics(alt, flatten=FALSE), perFeatureQCMetrics(original, flatten=FALSE))

    expect_equal(perCellQCMetrics(alt, subset=list(set=1:10), flatten=FALSE), 
        perCellQCMetrics(original, subset=list(set=1:10), flatten=FALSE))
    expect_equal(perFeatureQCMetrics(alt, subset=list(set=1:10), flatten=FALSE), 
        perFeatureQCMetrics(original, subset=list(set=1:10), flatten=FALSE))
})

test_that("we can compute standard QC metrics across multiple cores", {
    expect_equal(perCellQCMetrics(original, flatten=FALSE), 
        perCellQCMetrics(original, BPPARAM=safeBPParam(3), flatten=FALSE))
    expect_equal(perFeatureQCMetrics(original, flatten=FALSE), 
        perFeatureQCMetrics(original, BPPARAM=safeBPParam(3), flatten=FALSE))

    expect_equal(perCellQCMetrics(original, subset=list(set=1:10), flatten=FALSE), 
        perCellQCMetrics(original, subset=list(set=1:10), BPPARAM=safeBPParam(3), flatten=FALSE))
    expect_equal(perFeatureQCMetrics(original, subset=list(set=1:10), flatten=FALSE), 
        perFeatureQCMetrics(original, subset=list(set=1:10), BPPARAM=safeBPParam(3), flatten=FALSE))
})

#######################################################################
# addPerCellQCMetrics works: 

test_that("addPerCellQCMetrics adds rowData columns", {
    genes_numeric <- 1:5
    genes_names <- c("Gene_0001", "Gene_2000")
    genes_logical <- vector("logical", nrow(original))
    genes_logical[sample(nrow(original), 5)] <- TRUE
    
    subsets <- list(numeric = genes_numeric, logical = genes_logical, names = genes_names)
    out_sce <- addPerCellQCMetrics(original, subsets = subsets)
    rd <- rowData(out_sce)

    expect_identical(rd$subset_numeric, seq_len(nrow(original)) %in% subsets$numeric)
    expect_identical(rd$subset_logical, subsets$logical)
    expect_identical(rd$subset_names, rownames(original) %in% subsets$names)

    out_sce <- addPerCellQCMetrics(original, subsets = subsets, subset.prefix=NULL)
    rd <- rowData(out_sce)
    expect_null(rd$subset_numeric)
    expect_null(rd$subset_logical)
    expect_null(rd$subset_names)

    # Checking that the colData is modified
    expect_type(out_sce$sum, "double")
    expect_type(out_sce$detected, "double")
})
