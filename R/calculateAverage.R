#' Calculate per-feature average counts
#'
#' Calculate the average count for each feature after normalizing observations using per-cell size factors. 
#'
#' @param x A numeric matrix of counts where features are rows and columns are cells.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size.factors A numeric vector containing size factors.
#' If \code{NULL}, these are calculated or extracted from \code{x}.
#' @param assay.type A string specifying the assay of \code{x} containing the count matrix.
#' @param subset.row A vector specifying the subset of rows of \code{object} for which to return a result.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the calculations should be parallelized.
#' Only relevant for parallelized \code{\link{rowSums}(x)}, e.g., for \linkS4class{DelayedMatrix} inputs.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param size_factors,subset_row,exprs_values Soft-deprecated counterparts to the arguments above.
#'
#' @details 
#' The size factor-adjusted average count is defined by dividing each count by the size factor and taking the average across cells.
#' All size factors are scaled so that the mean is 1 across all cells, to ensure that the averages are interpretable on the same scale of the raw counts. 
#'
#' If no size factors are supplied, they are determined automatically:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size_factors} are supplied, they will override any size factors present in \code{x}.
#'
#' @return A numeric vector of average count values with same length as number of features 
#' (or the number of features in \code{subset_row} if supplied).
#' 
#' @author Aaron Lun
#'
#' @name calculateAverage
#'
#' @seealso
#' \code{\link{librarySizeFactors}}, for the default calculation of size factors.
#'
#' \code{\link{logNormCounts}}, for the calculation of normalized expression values.
#'
#' @examples
#' example_sce <- mockSCE()
#' ave_counts <- calculateAverage(example_sce)
#' summary(ave_counts)
NULL

#' @importFrom BiocParallel register bpparam SerialParam
#' @importFrom Matrix rowMeans
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.calculate_average <- function(x, size.factors=NULL, subset.row=NULL, BPPARAM = SerialParam(),
    size_factors=NULL, subset_row=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    size.factors <- .replace(size.factors, size_factors)

    subset.row <- .subset2index(subset.row, x, byrow=TRUE)
    size.factors <- .get_default_sizes(x, size.factors, center.size.factors=TRUE, subset.row=subset.row)

    # For DelayedArray's parallelized row/colSums.
    oldbp <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldbp))

    rowMeans(normalizeCounts(x, size.factors, subset.row=subset.row, log=FALSE))
}

#' @export
#' @rdname calculateAverage
setGeneric("calculateAverage", function(x, ...) standardGeneric("calculateAverage"))

#' @export
#' @rdname calculateAverage
setMethod("calculateAverage", "ANY", .calculate_average)

#' @export
#' @rdname calculateAverage
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateAverage", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) { 
    assay.type <- .replace(assay.type, exprs_values)
    .calculate_average(assay(x, assay.type), ...)
})

#' @export
#' @rdname calculateAverage
#' @importFrom BiocGenerics sizeFactors
#' @importFrom SingleCellExperiment altExp
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("calculateAverage", "SingleCellExperiment", function(x, size.factors=NULL, ...) { 
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    callNextMethod(x, size.factors=size.factors, ...)
})
