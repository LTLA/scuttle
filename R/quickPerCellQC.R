#' Quick cell-level QC
#'
#' A convenient utility that identifies low-quality cells based on frequently used QC metrics.
#'
#' @param x A \linkS4class{DataFrame} containing per-cell QC statistics, as computed by \code{\link{perCellQCMetrics}}.
#' Alternatively, a \linkS4class{SummarizedExperiment} object that can be used to create such a DataFrame via \code{\link{perCellQCMetrics}}.
#' @param sum.field String specifying the column of \code{x} containing the library size for each cell.
#' @param detected.field String specifying the column of \code{x} containing the number of detected features per cell.
#' @param sub.fields Character vector specifying the column(s) of \code{x} containing the percentage of counts in subsets of \dQuote{control features}, usually mitochondrial genes or spike-in transcripts.
#'
#' If set to \code{TRUE}, this will default to all columns in \code{x} with names following the patterns \code{"subsets_.*_percent"} and \code{"altexps_.*_percent"}. 
#' @param ... For the generic, further arguments to pass to specific methods.
#' 
#' For the ANY method, further arguments to pass to \code{\link{isOutlier}}.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param subsets,assay.type Arguments to pass to the \code{\link{perCellQCMetrics}} function, exposed here for convenience.
#' @param other.args A named list containing other arguments to pass to the \code{\link{perCellQCMetrics}} function.
#' @param lib_size,n_features,percent_subsets Soft-deprecated equivalents of the arguments above. 
#'
#' @return
#' A \linkS4class{DataFrame} with one row per cell and containing columns of logical vectors.
#' Each column specifies a reason for why a cell was considered to be low quality,
#' with the final \code{discard} column indicating whether the cell should be discarded.
#'
#' @details
#' The ANY method simply calls \code{\link{isOutlier}} on the various QC metrics in \code{x}.
#' \itemize{
#' \item For \code{sum.field}, small outliers are detected. 
#' These are considered to represent low-quality cells that have not been insufficiently sequenced.
#' Detection is performed on the log-scale to adjust for a heavy right tail and to improve resolution at zero.
#' \item For \code{detected.field}, small outliers are detected.
#' These are considered to represent low-quality cells with low-complexity libraries. 
#' Detection is performed on the log-scale to adjust for a heavy right tail.
#' This is done on the log-scale to adjust for a heavy right tail and to improve resolution at zero.
#' \item For each column specified by \code{sub.fields}, large outliers are detected.
#' This aims to remove cells with high spike-in or mitochondrial content, usually corresponding to damaged cells.
#' While these distributions often have heavy right tails, the putative low-quality cells are often present in this tail;
#' thus, transformation is not performed to ensure maintain resolution of the filter.
#' }
#'
#' Users can control the outlier detection (e.g., change the number of MADs, specify batches)
#' by passing appropriate arguments to \code{...}.
#'
#' The SummarizedExperiment method calls \code{\link{perCellQCMetrics}} before running outlier detection.
#' This is simply a convenient wrapper that avoids the need for two separate function calls in routine applications.
#'
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' x <- perCellQCMetrics(example_sce, subsets=list(Mito=1:100))
#'
#' discarded <- quickPerCellQC(x, percent_subsets=c(
#'     "subsets_Mito_percent", "altexps_Spikes_percent"))
#' colSums(as.data.frame(discarded))
#'
#' @seealso
#' \code{\link{perCellQCMetrics}}, for calculation of these metrics.
#'
#' \code{\link{isOutlier}}, to identify outliers with a MAD-based approach.
#' @export
#' @name quickPerCellQC
NULL

#' @importFrom S4Vectors DataFrame
.quickPerCellQC <- function(x, sum.field="sum", detected.field="detected", sub.fields=NULL, 
    ..., lib_size=NULL, n_features=NULL, percent_subsets=NULL)
{
    sum.field <- .replace(sum.field, lib_size)
    detected.field <- .replace(detected.field, n_features)
    sub.fields <- .replace(sub.fields, percent_subsets)

    output <- DataFrame(
        low_lib_size=isOutlier(x[[sum.field]], log=TRUE, type="lower", ...),
        low_n_features=isOutlier(x[[detected.field]], log=TRUE, type="lower", ...)
    )

    if (isTRUE(sub.fields)) {
        sub.fields <- colnames(x)[
            grepl("^subsets_.*_percent", colnames(x)) |
            grepl("^altexps_.*_percent", colnames(x))
        ]
    }

    for (i in sub.fields) {
        output[[paste0("high_", i)]] <- isOutlier(x[[i]], type="higher", ...)
    }

    discard <- Reduce("|", output)
    output$discard <- discard
    output
}

#' @export
#' @rdname quickPerCellQC
setGeneric("quickPerCellQC", function(x, ...) standardGeneric("quickPerCellQC"))

#' @export
#' @rdname quickPerCellQC
setMethod("quickPerCellQC", "ANY", .quickPerCellQC)

#' @export
#' @rdname quickPerCellQC
setMethod("quickPerCellQC", "SummarizedExperiment", function(x, ..., subsets=NULL, assay.type="counts", other.args=list()) {
    x <- do.call(perCellQCMetrics, c(list(x, subsets=subsets, assay.type=assay.type), other.args))
    .quickPerCellQC(x, ...)
})
