#' Quick cell-level QC
#'
#' A convenient utility that identifies low-quality cells based on frequently used QC metrics.
#'
#' @param df A \linkS4class{DataFrame} containing per-cell QC statistics, as computed by \code{\link{perCellQCMetrics}}.
#' @param sum.field String specifying the column of \code{df} containing the library size for each cell.
#' @param detected.field String specifying the column of \code{df} containing the number of detected features per cell.
#' @param subset.fields Character vector specifying the column(s) of \code{df} containing the percentage of counts in subsets of \dQuote{control features}, usually mitochondrial genes or spike-in transcripts.
#' @param ... Further arguments to pass to \code{\link{isOutlier}}.
#' @param lib_size,n_detected,percent_subset Soft-deprecated equivalents of the arguments above. 
#'
#' @return
#' A \linkS4class{DataFrame} with one row per cell and containing columns of logical vectors.
#' Each column specifies a reason for why a cell was considered to be low quality,
#' with the final \code{discard} column indicating whether the cell should be discarded.
#'
#' @details
#' This function simply calls \code{\link{isOutlier}} on the various QC metrics in \code{df}.
#' \itemize{
#' \item For \code{sum.field}, small outliers are detected. 
#' These are considered to represent low-quality cells that have not been insufficiently sequenced.
#' Detection is performed on the log-scale to adjust for a heavy right tail and to improve resolution at zero.
#' \item For \code{detected.field}, small outliers are detected.
#' These are considered to represent low-quality cells with low-complexity libraries. 
#' Detection is performed on the log-scale to adjust for a heavy right tail.
#' This is done on the log-scale to adjust for a heavy right tail and to improve resolution at zero.
#' \item For each column specified by \code{subset.fields}, large outliers are detected.
#' This aims to remove cells with high spike-in or mitochondrial content, usually corresponding to damaged cells.
#' While these distributions often have heavy right tails, the putative low-quality cells are often present in this tail;
#' thus, transformation is not performed to ensure maintain resolution of the filter.
#' }
#'
#' Users can control the outlier detection (e.g., change the number of MADs, specify batches)
#' by passing appropriate arguments to \code{...}.
#'
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' df <- perCellQCMetrics(example_sce, subsets=list(Mito=1:100))
#'
#' discarded <- quickPerCellQC(df, percent_subsets=c(
#'     "subsets_Mito_percent", "altexps_Spikes_percent"))
#' colSums(as.data.frame(discarded))
#'
#' @seealso
#' \code{\link{perCellQCMetrics}}, for calculation of these metrics.
#'
#' \code{\link{isOutlier}}, to identify outliers with a MAD-based approach.
#' @export
#' @importFrom S4Vectors DataFrame
quickPerCellQC <- function(df, sum.field="sum", detected.field="detected", subset.fields=NULL, 
    ..., lib_size=NULL, n_features=NULL, percent_subsets=NULL)
{
    sum.field <- .replace(sum.field, lib_size)
    detected.field <- .replace(detected.field, n_detected)
    subset.fields <- .replace(subset.fields, percent_subsets)

    output <- DataFrame(
        low_lib_size=isOutlier(df[[sum.field]], log=TRUE, type="lower", ...),
        low_n_features=isOutlier(df[[detected.field]], log=TRUE, type="lower", ...)
    )

    for (i in subset.fields) {
        output[[paste0("high_", i)]] <- isOutlier(df[[i]], type="higher", ...)
    }

    discard <- Reduce("|", output)
    output$discard <- discard
    output
}
