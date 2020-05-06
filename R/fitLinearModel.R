#' Fit a linear model
#'
#' No-frills fitting of a linear model to the rows of any matrix-like object with numeric values.
#' 
#' @param x A numeric matrix-like object where columns are samples (e.g., cells) and rows are usually features (e.g., genes).
#' @param design A numeric design matrix with number of rows equal to \code{ncol(x)}.
#' This should be of full column rank.
#' @param get.coefs A logical scalar indicating whether the coefficients should be returned.
#' @param subset.row An integer, character or logical vector indicating the rows of \code{x} to use for model fitting.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization backend to use.
#'
#' @return
#' If \code{get.coefs=TRUE}, a list is returned containing:
#' \itemize{
#' \item \code{coefficents}, a numeric matrix of coefficient estimates,
#' with one row per row of \code{x} (or a subset thereof specified by \code{subset.row}) 
#' and one column per column of \code{design}.
#' \item \code{mean}, a numeric vector of row means of \code{x}.
#' Computed as a courtesy to avoid iterating over the matrix twice.
#' \item \code{variance}, a numeric vector of residual variances per row of \code{x}.
#' Computed by summing the residual effects from the fitted model.
#' }
#' 
#' Otherwise, if \code{get.coefs=FALSE}, the same list is returned without \code{coefficients}.
#'
#' @details 
#' This function is basically a stripped-down version of \linkS4class{lm.fit},
#' made to operate on any matrix representation (ordinary, sparse, whatever).
#' It is generally intended for use inside other functions that require 
#' robust and efficient linear model fitting.
#'
#' @examples
#' y <- Matrix::rsparsematrix(1000, 1000, 0.1)
#' design <- model.matrix(~runif(1000))
#'
#' output <- fitLinearModel(y, design)
#' head(output$coefficients)
#' head(output$variance)
#'
#' @author Aaron Lun
#' @export
#' @importFrom BiocParallel bplapply SerialParam
fitLinearModel <- function(x, design, get.coefs=TRUE, subset.row=NULL, BPPARAM=SerialParam()) {
    QR <- .ranksafeQR(design)

    subset.row <- .subset2index(subset.row, x)
    by.rows <- .splitRowsByWorkers(x, BPPARAM=BPPARAM, subset.row=subset.row)
    bp.out <- bplapply(by.rows, FUN=fit_linear_model, qr=QR$qr, qraux=QR$qraux, get_coefs=get.coefs)

    all.means <- unlist(lapply(bp.out, "[[", i=2))
    all.vars <- unlist(lapply(bp.out, "[[", i=3))
    names(all.means) <- names(all.vars) <- rownames(x)
    output <- list(mean=all.means, variance=all.vars)

    if (get.coefs) {
        all.coefs <- do.call(cbind, lapply(bp.out, "[[", i=1))
        dimnames(all.coefs) <- list(colnames(design), rownames(x))
        output <- c(list(coefficients=t(all.coefs)), output)
    }

    output
}

#' @export
.ranksafeQR <- function(design, tol=1e-7) 
# Rank-checking QR decomposition of a design matrix. Throws an
# error if the design matrix is not of full rank, which simplifies
# downstream processes as full rank can always be assumed.
{
    out <- qr(design, LAPACK=TRUE)
    d <- diag(out$qr)
    if (!all(abs(d) > tol)) { 
        stop("design matrix is not of full rank")
    }
    out
}
