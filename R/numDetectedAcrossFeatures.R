#' Number of detected expression values per group of features
#' 
#' Computes the number of detected expression values (by default, defined as non-zero counts) 
#' for each group of features in each cell.
#'
#' @param x A numeric matrix of counts containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @inheritParams sumCountsAcrossFeatures
#' @param average Logical scalar indicating whether the proportion of non-zero counts in each group should be computed instead.
#' @param ... For the generic, further arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY method.
#' @param threshold A numeric scalar specifying the threshold above which a gene is considered to be detected.
#' 
#' @return An integer matrix containing the number of detected expression values in each group of features (row) and cell (column).
#' If \code{average=TRUE}, this is instead a numeric matrix containing the proportion of detected values.
#'
#' @author Aaron Lun
#' @seealso
#' \code{\link{sumCountsAcrossFeatures}}, on which this function is based.
#' 
#' @examples
#' example_sce <- mockSCE()
#'
#' ids <- sample(paste0("GENE_", 1:100), nrow(example_sce), replace=TRUE)
#' byrow <- numDetectedAcrossFeatures(example_sce, ids)
#' head(byrow[,1:10])
#'
#' @name numDetectedAcrossFeatures
NULL

#' @importFrom BiocParallel SerialParam 
.nexprs_across_features <- function(x, ids, subset.row=NULL, subset.col=NULL, 
    average=FALSE, threshold=0, BPPARAM=SerialParam(), 
    subset_row=NULL, subset_col=NULL, detection_limit=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    subset.col <- .replace(subset.col, subset_col)
    threshold <- .replace(threshold, detection_limit)

    .sum_across_features(x, ids, subset.row=subset.row, subset.col=subset.col, 
        average=average, BPPARAM=BPPARAM, modifier=function(x) x > threshold)
} 

#' @export
#' @rdname numDetectedAcrossFeatures
setGeneric("numDetectedAcrossFeatures", function(x, ...) standardGeneric("numDetectedAcrossFeatures"))

#' @export
#' @rdname numDetectedAcrossFeatures
setMethod("numDetectedAcrossFeatures", "ANY", .nexprs_across_features)

#' @export
#' @rdname numDetectedAcrossFeatures
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("numDetectedAcrossFeatures", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .nexprs_across_features(assay(x, assay.type), ...)
})
