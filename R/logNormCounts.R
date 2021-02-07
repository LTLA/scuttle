#' Compute log-normalized expression values
#'
#' Compute log-transformed normalized expression values from a count matrix in a \linkS4class{SingleCellExperiment} object.
#'
#' @param x A \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing a count matrix.
#' @inheritParams normalizeCounts
#' @param use.altexps,use_altexps Deprecated, use \code{\link{applySCE}} instead (see Examples).
#' @param ... For the generic, additional arguments passed to specific methods. 
#'
#' For the methods, additional arguments passed to \code{\link{normalizeCounts}}.
#' @param name String containing an assay name for storing the output normalized values.
#' Defaults to \code{"logcounts"} when \code{transform="log"}, \code{"ashcounts"} when \code{transform="asinh"}, and \code{"normcounts"} otherwise.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how library size factor calculations should be parallelized.
#' Only used if \code{size.factors} is not specified.
#'
#' @details
#' This function is a convenience wrapper around \code{\link{normalizeCounts}}.
#' It returns a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} containing the normalized values in a separate assay.
#' This makes it easier to perform normalization by avoiding book-keeping errors during a long analysis workflow.
#' 
#' If \code{NULL}, size factors are determined as described in \code{\link{normalizeCounts}}.
#' \code{subset.row} and \code{normalize.all} have the same interpretation as for \code{\link{normalizeCounts}}.
#'
#' @return 
#' \code{x} is returned containing the (log-)normalized expression values in an additional assay named as \code{name}.
#' 
#' If \code{x} is a \linkS4class{SingleCellExperiment}, the size factors used for normalization are stored in \code{\link{sizeFactors}}.
#' These are centered if \code{center.size.factors=TRUE}.
#'
#' @author Aaron Lun, based on code by Davis McCarthy 
#' @seealso
#' \code{\link{normalizeCounts}}, which is used to compute the normalized expression values.
#'
#' @examples
#' example_sce <- mockSCE()
#'
#' # Standard library size normalization:
#' example_sce2 <- logNormCounts(example_sce)
#' assayNames(example_sce2)
#' logcounts(example_sce2)[1:5,1:5]
#'
#' # Without logging, the assay is 'normcounts':
#' example_sce2 <- logNormCounts(example_sce, log=FALSE)
#' assayNames(example_sce2)
#' normcounts(example_sce2)[1:5,1:5]
#'
#' # Pre-loading with size factors:
#' example_sce2 <- computeMedianFactors(example_sce)
#' example_sce2 <- logNormCounts(example_sce2)
#' logcounts(example_sce2)[1:5,1:5]
#'
#' # Also normalizing the alternative experiments:
#' example_sce2 <- applySCE(example_sce, logNormCounts)
#' logcounts(altExp(example_sce2))[1:5,1:5]
#'
#' @name logNormCounts
NULL

#' @export
#' @rdname logNormCounts
setGeneric("logNormCounts", function(x, ...) standardGeneric("logNormCounts"))

#' @export
#' @rdname logNormCounts
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("logNormCounts", "SummarizedExperiment", function(x, size.factors=NULL, 
    log=NULL, transform=c("log", "none", "asinh"), pseudo.count=1, 
    center.size.factors=TRUE, ..., subset.row=NULL, normalize.all=FALSE, 
    assay.type="counts", name=NULL, BPPARAM=SerialParam(), 
    size_factors=NULL, pseudo_count=NULL, center_size_factors=NULL, exprs_values=NULL)
{
    size.factors <- .replace(size.factors, size_factors)
    pseudo.count <- .replace(pseudo.count, pseudo_count)
    center.size.factors <- .replace(center.size.factors, center_size_factors)
    assay.type <- .replace(assay.type, exprs_values)
    transform <- .choose_transform(log, match.arg(transform))

    if (!is.null(subset.row) && !normalize.all) {
        x <- x[subset.row,,drop=FALSE]
        subset.row <- NULL
    }
    
    FUN <- .se_lnc(assay.type=assay.type, transform=transform, pseudo.count=pseudo.count, ..., 
        normalize.all=normalize.all, subset.row=subset.row, name=name, BPPARAM=BPPARAM) 

    FUN(x, size.factors=size.factors, center.size.factors=center.size.factors)
})

#' @importFrom SummarizedExperiment assay<-
.se_lnc <- function(assay.type, transform, pseudo.count, ..., name) {
    args <- list(..., assay.type=assay.type, transform=transform, pseudo.count=pseudo.count)
    if (is.null(name)) {
        name <- switch(transform,
            none="normcounts",
            log="logcounts",
            asinh="ashcounts"
        )
    }

    function(x, ...) {
        out <- do.call(normalizeCounts, c(list(x, ...), args))
        assay(x, name) <- out
        x
    }
}

#' @export
#' @rdname logNormCounts
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom SingleCellExperiment altExp altExp<- int_metadata int_metadata<-
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("logNormCounts", "SingleCellExperiment", function(x, size.factors=sizeFactors(x), 
    log=NULL, transform=c("log", "none", "asinh"), pseudo.count=1, 
    center.size.factors=TRUE, ..., subset.row=NULL, normalize.all=FALSE, 
    assay.type="counts", use.altexps=FALSE, name=NULL, BPPARAM=SerialParam(), 
    size_factors=NULL, pseudo_count=NULL, center_size_factors=NULL, exprs_values=NULL, use_altexps=NULL) 
{
    size.factors <- .replace(size.factors, size_factors)
    pseudo.count <- .replace(pseudo.count, pseudo_count)
    center.size.factors <- .replace(center.size.factors, center_size_factors)
    use.altexps <- .replace(use.altexps, use_altexps)
    assay.type <- .replace(assay.type, exprs_values)
    transform <- .choose_transform(log, match.arg(transform))

    if (!is.null(subset.row) && !normalize.all) {
        x <- x[subset.row,,drop=FALSE]
        subset.row <- NULL
    }

    # Guarantee that we store (centered) size factors in the output matrix.
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, assay.type=assay.type, BPPARAM=BPPARAM, subset.row=subset.row)
    }
    size.factors <- .center.size.factors(size.factors, center.size.factors)
    sizeFactors(x) <- size.factors

    # Set center.size.factors=FALSE, as we've already centered above.
    FUN <- .se_lnc(assay.type=assay.type, transform=transform, pseudo.count=pseudo.count, ..., name=name) 
    x <- FUN(x, size.factors=size.factors, center.size.factors=FALSE)
    if (transform=="log") {
        if (is.null(int_metadata(x)$scater)) {
            int_metadata(x)$scater <- list()
        }
        int_metadata(x)$scater$pseudo.count <- pseudo.count
    }

    if (!is.null(use.altexps) && !isFALSE(use.altexps)) {
        .Deprecated(msg="'use.altexps=' is deprecated.\nUse 'applySCE(x, logNormCounts)' instead.")
        use.altexps <- .use_names_to_integer_indices(use.altexps, x=x, nameFUN=altExpNames, msg="use.altexps")
        for (i in use.altexps) {
            tryCatch({
                altExp(x, i) <- FUN(altExp(x, i), center.size.factors=center.size.factors)
            }, error=function(err) {
                stop(paste0(sprintf("failed to normalize 'altExp(x, %s)'\n", deparse(i)), conditionMessage(err)))
            })
        }
    }

    x
})
