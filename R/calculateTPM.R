#' Calculate TPMs
#'
#' Calculate transcripts-per-million (TPM) values for expression from feature-level counts.
#'
#' @param x A numeric matrix of counts where features are rows and cells are columns.
#'
#' Alternatively, a \linkS4class{SummarizedExperiment} or a \linkS4class{SingleCellExperiment} containing such counts.
#' @param size_factors A numeric vector containing size factors to adjust the library sizes.
#' If \code{NULL}, the library sizes are used directly. 
#' @param lengths Numeric vector providing the effective length for each feature in \code{x}.
#' Alternatively \code{NULL}, see Details.
#' @param exprs_values A string specifying the assay of \code{x} containing the count matrix.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the ANY method, further arguments to pass to \code{\link{calculateCPM}}.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#'
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#'
#' @details
#' For read count data, this function assumes uniform coverage along the (effective) length of the transcript.
#' Thus, the number of transcripts for a gene is proportional to the read count divided by the transcript length.
#' Here, the division is done before calculation of the library size to compute per-million values,
#' where \code{\link{calculateFPKM}} will only divide by the length after library size normalization.
#'
#' For UMI count data, this function should be run with \code{lengths=NULL}, i.e., no division by the effective length.
#' This is because the number of UMIs is a direct (albeit biased) estimate of the number of transcripts.
#'
#' @return A numeric matrix of TPM values with the same dimensions as \code{x} (unless \code{subset.row} is specified).
#'
#' @author Aaron Lun, based on code by Davis McCarthy
#' @seealso
#' \code{\link{calculateCPM}}, on which this function is based.
#'
#' @examples
#' example_sce <- mockSCE()
#' eff_len <- runif(nrow(example_sce), 500, 2000)
#' tout <- calculateTPM(example_sce, lengths = eff_len)
#' str(tout)
#'
#' @name calculateTPM
NULL

.calculate_tpm <- function(x, lengths=NULL, ...) {
    if (!is.null(lengths)) {
        x <- x/lengths
    }
    .calculate_cpm(x, ...)
}

#' @export
#' @rdname calculateTPM
setGeneric("calculateTPM", function(x, ...) standardGeneric("calculateTPM"))

#' @export
#' @rdname calculateTPM
setMethod("calculateTPM", "ANY", .calculate_tpm)

#' @export
#' @rdname calculateTPM
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("calculateTPM", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .calculate_tpm(assay(x, assay.type), ...)
})

#' @export
#' @rdname calculateTPM
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment 
setMethod("calculateTPM", "SingleCellExperiment", function(x, lengths=NULL, size.factors=NULL, ...) {
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    callNextMethod(x=x, lengths=lengths, size.factors=size.factors, ...)
})
