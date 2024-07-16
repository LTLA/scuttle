#' Compute library size factors
#' 
#' Define per-cell size factors from the library sizes (i.e., total sum of counts per cell).
#'
#' @param x For \code{librarySizeFactors}, a numeric matrix of counts with one row per feature and column per cell.
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} containing such counts.
#'
#' For \code{computeLibraryFactors}, only a \linkS4class{SingleCellExperiment} containing a count matrix is accepted.
#' @param subset.row A vector specifying whether the size factors should be computed from a subset of rows of \code{x}.
#' @param assay.type String or integer scalar indicating the assay of \code{x} containing the counts.
#' @param geometric Deprecated, logical scalar indicating whether the size factor should be defined using the geometric mean.
#' @param pseudo_count Deprecated, numeric scalar specifying the pseudo-count to add when \code{geometric=TRUE}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how calculations are to be parallelized.
#' Only relevant when \code{x} is a \linkS4class{DelayedArray} object.
#' @param ... For the \code{librarySizeFactors} generic, arguments to pass to specific methods.
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For \code{computeLibraryFactors}, further arguments to pass to \code{librarySizeFactors}.
#' @param subset_row,exprs_values Soft-deprecated equivalents to the arguments above.
#'
#' @details
#' Library sizes are converted into size factors by scaling them so that their mean across cells is unity.
#' This ensures that the normalized values are still on the same scale as the raw counts.
#' Preserving the scale is useful for interpretation of operations on the normalized values,
#' e.g., the pseudo-count used in \code{\link{logNormCounts}} can actually be considered an additional read/UMI.
#' This is important for ensuring that the effect of the pseudo-count decreases with increasing sequencing depth,
#' see \code{?\link{normalizeCounts}} for a discussion of this effect.
#'
#' With library size-derived size factors, we implicitly assume that sequencing coverage is the only difference between cells.
#' This is reasonable for homogeneous cell populations but is compromised by composition biases from DE between cell types.
#' In such cases, the library size factors will not be correct though any effects on downstream conclusions will vary,
#' e.g., clustering is usually unaffected by composition biases but log-fold change estimates will be less accurate.
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
#' For \code{librarySizeFactors}, a numeric vector of size factors is returned for all methods.
#'
#' For \code{computeLibraryFactors}, \code{x} is returned containing the size factors in \code{\link{sizeFactors}(x)}.
#'
#' @name librarySizeFactors
#' @examples
#' example_sce <- mockSCE()
#' summary(librarySizeFactors(example_sce))
NULL

#' @importFrom BiocParallel SerialParam
#' @importFrom MatrixGenerics colSums colMeans
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.library_size_factors <- function(x, subset.row=NULL, geometric=FALSE, BPPARAM=SerialParam(),
    subset_row=NULL, pseudo_count=1) 
{
    subset.row <- .replace(subset.row, subset_row)

    if (!is.null(subset.row)) {
       x <- x[subset.row,,drop=FALSE]
    }

    oldBP <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldBP))

    if (!geometric) {
        lib.sizes <- colSums(x)
        lib.sizes/mean(lib.sizes)
    } else {
        .Deprecated(new="geometricSizeFactors")
        geometricSizeFactors(x, subset.row=subset.row, pseudo.count=pseudo_count)
    }
}       

#' @export
#' @rdname librarySizeFactors
setGeneric("librarySizeFactors", function(x, ...) standardGeneric("librarySizeFactors"))

#' @export
#' @rdname librarySizeFactors
setMethod("librarySizeFactors", "ANY", .library_size_factors)

#' @export
#' @rdname librarySizeFactors
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("librarySizeFactors", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .library_size_factors(assay(x, assay.type), ...)
})

#' @export
#' @rdname librarySizeFactors
#' @importFrom BiocGenerics sizeFactors<-
computeLibraryFactors <- function(x, ...) {
    sf <- librarySizeFactors(x, ...)
    sizeFactors(x) <- sf
    x
}
