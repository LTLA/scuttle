#' Per-feature quality control metrics
#'
#' Compute per-feature quality control metrics for a count matrix or a \linkS4class{SummarizedExperiment}.
#'
#' @param subsets A named list containing one or more vectors 
#' (a character vector of cell names, a logical vector, or a numeric vector of indices),
#' used to identify interesting sample subsets such as negative control wells.
#' @inheritParams perCellQCMetrics
#' @param detection_limit,exprs_values Soft deprecated equivalents to the arguments described above.
#'
#' @return
#' A \linkS4class{DataFrame} of QC statistics where each row corresponds to a row in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{mean}: numeric, the mean counts for each feature.
#' \item \code{detected}: numeric, the percentage of observations above \code{threshold}.
#' }
#'
#' If \code{flatten=FALSE}, the output DataFrame also contains the \code{subsets} field.
#' This a nested DataFrame containing per-feature QC statistics for each subset of columns.
#'
#' If \code{flatten=TRUE}, \code{subsets} is flattened to remove the hierarchical structure.
#' 
#' @author Aaron Lun
#' 
#' @details
#' This function calculates useful QC metrics for features, including the mean across all cells
#' and the number of expressed features (i.e., counts above the detection limit).
#' 
#' If \code{subsets} is specified, the same statistics are computed for each subset of cells.
#' This is useful for obtaining statistics for cell sets of interest, e.g., negative control wells.
#' These statistics are stored as nested \linkS4class{DataFrame}s in the output.
#' For example, if \code{subsets} contained \code{"empty"} and \code{"cellpool"}, the output would look like:
#' \preformatted{  output 
#'   |-- mean 
#'   |-- detected
#'   +-- subsets
#'       |-- empty
#'       |   |-- mean 
#'       |   |-- detected
#'       |   +-- ratio
#'       +-- cellpool 
#'           |-- mean
#'           |-- detected
#'           +-- ratio
#' }
#' The \code{ratio} field contains the ratio of the mean within each subset to the mean across all cells.
#' 
#' If \code{flatten=TRUE}, the nested DataFrames are flattened by concatenating the column names with underscores.
#' This means that, say, the \code{subsets$empty$mean} nested field becomes the top-level \code{subsets_empty_mean} field.
#' A flattened structure is more convenient for end-users performing interactive analyses,
#' but less convenient for programmatic access as artificial construction of strings is required.
#' @examples
#' example_sce <- mockSCE()
#' stats <- perFeatureQCMetrics(example_sce)
#' stats
#'
#' # With subsets.
#' stats2 <- perFeatureQCMetrics(example_sce, subsets=list(Empty=1:10))
#' stats2
#'
#' @seealso 
#' \code{\link{addPerFeatureQCMetrics}}, to add the QC metrics to the row metadata.
#' @export
#' @name perFeatureQCMetrics
NULL

#' @importFrom beachmat rowBlockApply
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom BiocParallel bplapply SerialParam
#' @importClassesFrom S4Vectors DFrame 
.per_feature_qc_metrics <- function(x, subsets = NULL, threshold = 0, BPPARAM=SerialParam(), flatten=TRUE,
    detection_limit=NULL) 
{
    threshold <- .replace(threshold, detection_limit)

    if (length(subsets) && is.null(names(subsets))){ 
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index, target = x, byrow = FALSE)

    # Computing all QC metrics, with cells split across workers.
    bp.out <- rowBlockApply(x, FUN=.per_feature_qc, cellcon=subsets, limit=threshold, BPPARAM=BPPARAM)

    # Aggregating across cores.
    full.info <- DataFrame(
        mean=unlist(lapply(bp.out, FUN=function(x) x[[1]][[1]])),
        detected=unlist(lapply(bp.out, FUN=function(x) x[[1]][[2]])) * 100,
        row.names=rownames(x)
    )

    # Collecting subset information.
    if (!is.null(subsets)) {
        sub.info <- make_zero_col_DFrame(nrow(x))
        for (i in seq_along(subsets)) {
            sub.out <- DataFrame(
                mean=unlist(lapply(bp.out, FUN=function(x) x[[2]][[i]][[1]])),
                detected=unlist(lapply(bp.out, FUN=function(x) x[[2]][[i]][[2]])) * 100
            )
            sub.out$ratio <- sub.out$mean/full.info$mean
            sub.info[[names(subsets)[i]]] <- sub.out
        }
        full.info$subsets <- sub.info
    }

    if (flatten) {
        full.info <- .flatten_nested_dims(full.info)
    }
    full.info
}

#' @importFrom MatrixGenerics rowMeans
#' @importClassesFrom SparseArray COO_SparseMatrix SVT_SparseMatrix
.per_feature_qc <- function(x, cellcon, limit) {
    if (is(x, "COO_SparseMatrix")) {
        x <- as(x, "SVT_SparseMatrix")
    }

    detected <- x > limit

    full <- list(
        sum=unname(rowMeans(x)),
        detected=unname(rowMeans(detected))
    )

    cellcons <- lapply(cellcon, function(i) {
        list(
            sum=unname(rowMeans(x[,i,drop=FALSE])),
            detected=unname(rowMeans(detected[,i,drop=FALSE]))
        )
    })

    list(full, cellcons)
}

#' @export
#' @rdname perFeatureQCMetrics
setGeneric("perFeatureQCMetrics", function(x, ...) standardGeneric("perFeatureQCMetrics"))

#' @export
#' @rdname perFeatureQCMetrics
setMethod("perFeatureQCMetrics", "ANY", .per_feature_qc_metrics)

#' @export
#' @rdname perFeatureQCMetrics
#' @importFrom SummarizedExperiment assay
setMethod("perFeatureQCMetrics", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .per_feature_qc_metrics(assay(x, assay.type), ...)
})
