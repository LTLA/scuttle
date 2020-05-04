#' Number of detected expression values per group of cells
#' 
#' Computes the number of detected expression values (by default, defined as non-zero counts) 
#' for each feature in each group of cells.
#'
#' @param x A numeric matrix of counts containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @inheritParams sumCountsAcrossCells
#' @param average Logical scalar indicating whether the proportion of non-zero counts in each group should be computed instead.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param threshold A numeric scalar specifying the threshold above which a gene is considered to be detected.
#' @param exprs_values,detection_limit Soft-deprecated equivalents of the arguments above.
#' 
#' @return 
#' A SummarizedExperiment is returned containing a count matrix in the first assay.
#' Each column corresponds to group as defined by a unique level or combination of levels in \code{ids}.
#' Each entry of the matrix contains the number of cells with detected expression for a feature and group.
#' The identities of the levels for each column are reported in the \code{\link{colData}}.
#' If \code{average=TRUE}, the assay is instead a numeric matrix containing the proportion of detected values.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{sumCountsAcrossCells}}, which computes the sum of counts within a group.
#' 
#' @examples
#' example_sce <- mockSCE()
#'
#' ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)
#' bycol <- numDetectedAcrossCells(example_sce, ids)
#' head(bycol)
#'
#' @name numDetectedAcrossCells
NULL

#' @importFrom BiocParallel SerialParam 
#' @importClassesFrom BiocParallel MulticoreParam
.nexprs_across_cells <- function(x, ids, subset.row=NULL, subset.col=NULL, 
    store.number="ncells", average=FALSE, threshold=0, BPPARAM=SerialParam(),
    subset_row=NULL, subset_col=NULL, store_number=NULL, detection_limit=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    subset.col <- .replace(subset.col, subset_col)
    store.number <- .replace(store.number, store_number)
    threshold <- .replace(threshold, detection_limit)

    .sum_counts_across_cells(x=x, ids=ids,subset.row=subset.row, subset.col=subset.col, 
        average=average, store.number=store.number, BPPARAM=BPPARAM, 
        modifier=function(x) (x > threshold) + 0L) # coercing to numeric to make life easier.
} 

#' @export
#' @rdname numDetectedAcrossCells
setGeneric("numDetectedAcrossCells", function(x, ...) standardGeneric("numDetectedAcrossCells"))

#' @export
#' @rdname numDetectedAcrossCells
setMethod("numDetectedAcrossCells", "ANY", .nexprs_across_cells)

#' @export
#' @rdname numDetectedAcrossCells
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("numDetectedAcrossCells", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .nexprs_across_cells(assay(x, assay.type), ...)
})
