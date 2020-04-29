#' Calculate CPMs
#'
#' Calculate counts-per-million (CPM) values from the count data.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size.factors A numeric vector containing size factors to adjust the library sizes.
#' If \code{NULL}, the library sizes are used directly. 
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the count matrix.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param size_factors,subset_row,exprs_values Soft-deprecated counterparts to the arguments above.
#'
#' @details 
#' If \code{size.factors} are provided or available in \code{x}, they are used to define the effective library sizes. 
#' This is done by scaling all size factors such that the mean factor is equal to the mean sum of counts across all features. 
#' The effective library sizes are then used as the denominator of the CPM calculation.
#'
#' @return A numeric matrix of CPM values with the same dimensions as \code{x} (unless \code{subset.row} is specified).
#'
#' @name calculateCPM
#' @author Aaron Lun
#' @seealso 
#' \code{\link{normalizeCounts}}, on which this function is based.
#'
#' @examples
#' example_sce <- mockSCE()
#' cpm(example_sce) <- calculateCPM(example_sce)
#' str(cpm(example_sce))
NULL

#' @importFrom Matrix colSums
.calculate_cpm <- function(x, size.factors=NULL, subset.row=NULL, size_factors=NULL, subset_row=NULL) {
    size.factors <- .replace(size.factors, size_factors)
    subset.row <- .replace(subset.row, subset_row)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    lib.sizes <- colSums(x) / 1e6
    if (!is.null(size.factors)) {
        lib.sizes <- size.factors / mean(size.factors) * mean(lib.sizes)
    }

    normalizeCounts(x, size.factors=lib.sizes, log=FALSE, center.size.factors=FALSE)
}

#' @export
#' @rdname calculateCPM
setGeneric("calculateCPM", function(x, ...) standardGeneric("calculateCPM"))

#' @export
#' @rdname calculateCPM
setMethod("calculateCPM", "ANY", .calculate_cpm)

#' @export
#' @rdname calculateCPM
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateCPM", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .calculate_cpm(assay(x, assay.type), ...)
})

#' @export
#' @rdname calculateCPM
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("calculateCPM", "SingleCellExperiment", function(x, size.factors=NULL, ...) {
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    callNextMethod(x=x, size.factors=size.factors, ...)
})
