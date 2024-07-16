#' Compute geometric size factors
#' 
#' Define per-cell size factors from the geometric mean of counts per cell.
#'
#' @param x For \code{geometricSizeFactors}, a numeric matrix of counts with one row per feature and column per cell.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such counts.
#'
#' For \code{computeGeometricFactors}, only a \linkS4class{SingleCellExperiment} containing a count matrix is accepted.
#' @param subset.row A vector specifying whether the size factors should be computed from a subset of rows of \code{x}.
#' @param assay.type String or integer scalar indicating the assay of \code{x} containing the counts.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add during log-transformation. 
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how calculations are to be parallelized.
#' Only relevant when \code{x} is a \linkS4class{DelayedArray} object.
#' @param ... For the \code{geometricSizeFactors} generic, arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{computeGeometricFactors}, further arguments to pass to \code{geometricSizeFactors}.
#'
#' @details
#' The geometric mean provides an alternative measure of the average coverage per cell,
#' in contrast to the library size factors (i.e., the arithmetic mean) computed by \code{\link{librarySizeFactors}}.
#' The main advantage of the geometric mean is that it is more robust to large outliers, due to the slowly increasing nature of the log-transform at large values;
#' in the normalization context, this translates to greater resistance to coposition biases from a few strongly upregulated genes.
#'
#' On the other hand, the geometric mean is a poor estimator of the relative bias at low or zero counts.
#' This is because the scaling of the coverage applies to the expectation of the raw counts, 
#' so the geometric mean only becomes an accurate estimator if the mean of the logs approaches the log of the mean (usually at high counts).
#' The arbitrary pseudo-count also has a bigger influence at low counts.
#'
#' As such, the geometric mean is only well-suited for deeply sequenced features, e.g., antibody-derived tags.
#' 
#' @author Aaron Lun
#'
#' @seealso 
#' \code{\link{normalizeCounts}} and \code{\link{logNormCounts}}, where these size factors are used by default.
#'
#' \code{\link{geometricSizeFactors}} and \code{\link{medianSizeFactors}}, 
#' for two other simple methods of computing size factors.
#' 
#' @return 
#' For \code{geometricSizeFactors}, a numeric vector of size factors is returned for all methods.
#'
#' For \code{computeGeometricFactors}, \code{x} is returned containing the size factors in \code{\link{sizeFactors}(x)}.
#'
#' @name geometricSizeFactors
#' @examples
#' example_sce <- mockSCE()
#' summary(geometricSizeFactors(example_sce))
NULL

#' @importFrom BiocParallel SerialParam
#' @importFrom MatrixGenerics colMeans
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.geometric_size_factors <- function(x, subset.row=NULL, pseudo.count=1, BPPARAM=SerialParam()) {
    if (!is.null(subset.row)) {
       x <- x[subset.row,,drop=FALSE]
    }

    oldBP <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldBP))

    geo <- 2^colMeans(normalizeCounts(x, size_factors=rep(1, ncol(x)), log=TRUE, pseudo.count=pseudo.count))
    geo/mean(geo)
}

#' @export
#' @rdname geometricSizeFactors
setGeneric("geometricSizeFactors", function(x, ...) standardGeneric("geometricSizeFactors"))

#' @export
#' @rdname geometricSizeFactors
setMethod("geometricSizeFactors", "ANY", .geometric_size_factors)

#' @export
#' @rdname geometricSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("geometricSizeFactors", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .geometric_size_factors(assay(x, assay.type), ...)
})

#' @export
#' @rdname geometricSizeFactors
#' @importFrom BiocGenerics sizeFactors<-
computeGeometricFactors <- function(x, ...) {
    sf <- geometricSizeFactors(x, ...)
    sizeFactors(x) <- sf
    x
}
