#' Quick cell-level QC
#'
#' A convenient utility that identifies low-quality cells based on frequently used QC metrics.
#'
#' @param x A \linkS4class{DataFrame} containing per-cell QC statistics, as computed by \code{\link{perCellQCMetrics}}.
#' Alternatively, a \linkS4class{SummarizedExperiment} object that can be used to create such a DataFrame via \code{\link{perCellQCMetrics}}.
#' @inheritParams perCellQCFilters 
#' @param ... For the generic, further arguments to pass to specific methods.
#' 
#' For the ANY method, further arguments to pass to \code{\link{isOutlier}}.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param subsets,assay.type Arguments to pass to the \code{\link{perCellQCMetrics}} function, exposed here for convenience.
#' @param filter Logical scalar indicating whether to filter out low-quality cells from \code{x}.
#' @param other.args A named list containing other arguments to pass to the \code{\link{perCellQCMetrics}} function.
#' @param lib_size,n_features,percent_subsets Soft-deprecated equivalents of the arguments above. 
#'
#' @return
#' If \code{filter=FALSE} or \code{x} is a DataFrame, a \linkS4class{DataFrame} is returned with one row per cell and containing columns of logical vectors.
#' Each column specifies a reason for why a cell was considered to be low quality,
#' with the final \code{discard} column indicating whether the cell should be discarded.
#'
#' If \code{filter=TRUE}, \code{x} is returned with the low-quality cells removed.
#' QC statistics and filtering information for all remaining cells are stored in the \code{\link{colData}}.
#'
#' @details
#' For DataFrame \code{x}, this function simply calls \code{\link{perCellQCFilters}}.
#' The latter should be directly used in such cases; DataFrame inputs are soft-deprecated here.
#'
#' For SummarizedExperiment \code{x}, this function is simply a convenient wrapper around \code{\link{perCellQCMetrics}} and \code{\link{perCellQCFilters}}.
#'
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#'
#' filtered_sce <- quickPerCellQC(example_sce, subsets=list(Mito=1:100),
#'     sub.fields=c("subsets_Mito_percent", "altexps_Spikes_percent"))
#' ncol(filtered_sce)
#'
#' # Same result as the longer chain of expressions:
#' stats <- perCellQCMetrics(example_sce, subsets=list(Mito=1:100))
#' discard <- perCellQCFilters(stats, 
#'     sub.fields=c("subsets_Mito_percent", "altexps_Spikes_percent"))
#' filtered_sce2 <- example_sce[,!discard$discard]
#' ncol(filtered_sce2)
#'
#' @seealso
#' \code{\link{perCellQCMetrics}}, for calculation of these metrics.
#'
#' \code{\link{perCellQCFilters}}, to define filter thresholds based on those metrics.
#' @export
#' @name quickPerCellQC
NULL

#' @export
#' @rdname quickPerCellQC
setGeneric("quickPerCellQC", function(x, ...) standardGeneric("quickPerCellQC"))

#' @export
#' @rdname quickPerCellQC
setMethod("quickPerCellQC", "ANY", function(x, sum.field="sum", detected.field="detected", sub.fields=NULL, 
    ..., lib_size=NULL, n_features=NULL, percent_subsets=NULL) 
{
    sum.field <- .replace(sum.field, lib_size)
    detected.field <- .replace(detected.field, n_features)
    sub.fields <- .replace(sub.fields, percent_subsets)
    perCellQCFilters(x, sum.field=sum.field, detected.field=detected.field, sub.fields=sub.fields, ...)
})

#' @export
#' @rdname quickPerCellQC
setMethod("quickPerCellQC", "SummarizedExperiment", function(x, ..., subsets=NULL, assay.type="counts", other.args=list(), filter=TRUE) {
    stats <- do.call(perCellQCMetrics, c(list(x, subsets=subsets, assay.type=assay.type), other.args))
    keep <- quickPerCellQC(stats, ...)
    if (!filter) {
        keep
    } else {
        colData(x) <- cbind(colData(x), stats, keep)
        x[,!keep$discard]
    }
})
