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
#' @param rank.error Logical scalar indicating whether to throw an error when \code{design} is not of full rank.
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
#' \item \code{residual.df}, an integer scalar containing the residual degrees of freedom for \code{design}.
#' }
#' 
#' Otherwise, if \code{get.coefs=FALSE}, the same list is returned without \code{coefficients}.
#'
#' @details 
#' This function is basically a stripped-down version of \code{\link{lm.fit}},
#' made to operate on any matrix representation (ordinary, sparse, whatever).
#' It is generally intended for use inside other functions that require robust and efficient linear model fitting.
#' 
#' If \code{design} is not of full rank and \code{rank.error=TRUE}, an error is raised.
#' If \code{rank.error=FALSE}, \code{NA} values are returned for all entries of the output list.
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
#' @importFrom beachmat rowBlockApply
fitLinearModel <- function(x, design, get.coefs=TRUE, subset.row=NULL, BPPARAM=SerialParam(), rank.error=TRUE) {
    QR <- .ranksafeQR(design, error=rank.error)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    if (is.null(QR)) {
        N <- nrow(x)
        dummy <- rep(NA_real_, N)
        names(dummy) <- rownames(x)
        output <- list(mean=dummy, variance=dummy, residual.df=NA_integer_)

        if (get.coefs) {
            mat <- matrix(NA_real_, N, ncol(design), dimnames=list(rownames(x), colnames(design)))
            output <- c(list(coefficients=mat), output)
        }

    } else {
        bp.out <- rowBlockApply(x, FUN=fit_linear_model, qr=QR$qr, qraux=QR$qraux, get_coefs=get.coefs, BPPARAM=BPPARAM)

        all.means <- unlist(lapply(bp.out, "[[", i=2))
        all.vars <- unlist(lapply(bp.out, "[[", i=3))
        names(all.means) <- names(all.vars) <- rownames(x)

        resid.df <- nrow(design) - ncol(design) # guaranteed to be full rank at this point.
        output <- list(mean=all.means, variance=all.vars, residual.df=resid.df)

        if (get.coefs) {
            all.coefs <- do.call(cbind, lapply(bp.out, "[[", i=1))
            all.coefs[QR$pivot,] <- all.coefs
            dimnames(all.coefs) <- list(colnames(design), rownames(x))
            output <- c(list(coefficients=t(all.coefs)), output)
        }
    }

    output
}

#' @export
.ranksafeQR <- function(design, tol=1e-7, error=TRUE) 
# Rank-checking QR decomposition of a design matrix. Throws an
# error if the design matrix is not of full rank, which simplifies
# downstream processes as full rank can always be assumed.
{
    out <- qr(design, LAPACK=TRUE)
    d <- diag(out$qr)
    if (!all(abs(d) > tol)) { 
        if (error) {
            stop("design matrix is not of full rank")
        } else {
            NULL
        }
    } else {
        out
    }
}
