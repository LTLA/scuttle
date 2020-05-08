# This tests the fitLinearModel function.
# library(scuttle); library(testthat); source("test-linear-model.R")

y <- matrix(rnorm(10000), ncol=20)
g <- gl(4, 5)
g2 <- rep(gl(2, 2), 5)
u <- runif(ncol(y))

test_that("fitLinearModel works on a variety of design matrices", {
    for (design in list(
            cbind(BLAH=rep(1, ncol(y))),
            model.matrix(~g),
            model.matrix(~g2),
            model.matrix(~u),
            model.matrix(~g + g2),
            model.matrix(~g + u),
            model.matrix(~g2 + u),
            model.matrix(~g2 + g),
            model.matrix(~u + g),
            model.matrix(~u + g2),
            model.matrix(~g + g2 + u),
            model.matrix(~u + g + g2),
            model.matrix(~g2 + u + g)
        )
    ) {
        fit <- fitLinearModel(y, design)
        ref <- lm.fit(x=design, t(y))

        expect_equal(fit$mean, rowMeans(y))
        expect_equal(fit$coefficients, t(ref$coefficients))
        expect_equal(fit$variance, colMeans(ref$effects[-seq_len(ref$rank),,drop=FALSE]^2))
    }
})

test_that("fitLinearModel works with the options", {
    design <- model.matrix(~g + g2 + u)
    ref <- fitLinearModel(y, design)

    sub <- fitLinearModel(y, design, subset.row=50:1)
    expect_identical(sub$coefficients, ref$coefficients[50:1,,drop=FALSE])
    expect_identical(sub$mean, ref$mean[50:1])
    expect_identical(sub$variance, ref$variance[50:1])

    par <- fitLinearModel(y, design, BPPARAM=BiocParallel::SnowParam(3))
    expect_identical(ref, par)
    par <- fitLinearModel(y, design, BPPARAM=safeBPParam(2))
    expect_identical(ref, par)

    nocoef <- fitLinearModel(y, design, get.coefs=FALSE)
    expect_identical(nocoef$mean, ref$mean)
    expect_identical(nocoef$variance, ref$variance)
    expect_identical(nocoef$coefficients, NULL)
})

test_that("fitLinearModel works with names", {
    design <- model.matrix(~g + g2 + u)
    rownames(y) <- paste0("GENE_", seq_len(nrow(y)))
    ref <- fitLinearModel(y, design)

    expect_identical(names(ref$mean), rownames(y))
    expect_identical(names(ref$variance), rownames(y))
    expect_identical(rownames(ref$coefficient), rownames(y))

    # Handles subsetting.
    chosen <- rbinom(nrow(y), 1, 0.5)==1
    sub <- fitLinearModel(y, design, subset.row=chosen)

    expect_identical(names(sub$mean), rownames(y)[chosen])
    expect_identical(names(sub$variance), rownames(y)[chosen])
    expect_identical(rownames(sub$coefficient), rownames(y)[chosen])
})

test_that("fitLinearModel correctly identifies low-rank matrices", {
    design <- model.matrix(~g) 
    design2 <- cbind(design, design)
    expect_error(fitLinearModel(y, design2), "not of full rank")
})

