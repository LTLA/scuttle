#' Calculate FPKMs
#'
#' Calculate fragments per kilobase of exon per million reads mapped (FPKM) values from the feature-level counts.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param lengths Numeric vector providing the effective length for each feature in \code{x}.
#' @param ... Further arguments to pass to \code{\link{calculateCPM}}.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param subset_row Soft-deprecated equivalent to the argument above.
#'
#' @details
#' FPKMs are computed by dividing the CPMs by the effective length of each gene in kilobases.
#' For RNA-seq datasets, the effective length is best set to the sum of lengths of all exons;
#' for nucleus sequencing datasets, the effective length may instead be the entire width of the gene body.
#' 
#' @return A numeric matrix of FPKM values with the same dimensions as \code{x} (unless \code{subset.row} is specified).
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#'
#' @seealso 
#' \code{\link{calculateCPM}}, for the initial calculation of CPM values.
#'
#' @examples
#' example_sce <- mockSCE()
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' fout <- calculateFPKM(example_sce, eff_len)
#' str(fout)
#' @export
calculateFPKM <- function(x, lengths, ..., subset.row=NULL, subset_row=NULL) {
    subset.row <- .replace(subset.row, subset_row)

    if (!is.null(subset.row)) {
        subset.row <- .subset2index(subset.row, x, byrow=TRUE)
        lengths <- lengths[subset.row]
    }

    out <- calculateCPM(x, subset.row=subset.row, ...)
    out / (lengths / 1e3)
}
