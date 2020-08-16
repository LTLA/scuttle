# This tests the quickPerCellQC function.
# library(testthat); library(scuttle); source("setup.R"); source("test-qc-quick.R")

original <- sce

test_that("quickPerCellQC works correctly", {
    df <- perCellQCMetrics(original)
    df$sum[1] <- 0
    df$detected[2] <- 0

    out <- quickPerCellQC(df)
    expect_true(out$low_lib_size[1])
    expect_true(out$low_n_features[2])

    expect_identical(out$discard, out$low_lib_size | out$low_n_features)
})

test_that("quickPerCellQC works correctly with subsets", {
    df <- perCellQCMetrics(original)
    df$subsets_BLAH_percent <- c(1, rep(0, nrow(df)-1))
    df$altexps_WHEE_percent <- c(0, 1, rep(0, nrow(df)-2))
    df$sum[3] <- 0
    df$detected[4] <- 0

    out <- quickPerCellQC(df, sub.fields=c("subsets_BLAH_percent", "altexps_WHEE_percent"))
    expect_true(out$high_subsets_BLAH_percent[1])
    expect_true(out$high_altexps_WHEE_percent[2])
    expect_identical(out$discard, Reduce("|", out[,1:4]))

    out2 <- quickPerCellQC(df, sub.fields=TRUE)
    expect_identical(out, out2)
})

test_that("quickPerCellQC works correctly with SCEs", {
    df <- perCellQCMetrics(original)
    out <- quickPerCellQC(df)
    expect_identical(out, quickPerCellQC(original))
})
