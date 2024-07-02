#' Add QC metrics to a SummarizedExperiment
#'
#' Convenient utilities to compute QC metrics and add them to a \linkS4class{SummarizedExperiment}'s row or column metadata.
#'
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
#' @inheritParams perFeatureQCMetrics
#' @param ... For \code{addPerCellQCMetrics}, further arguments to pass to \code{\link{perCellQCMetrics}}.
#' 
#' For \code{addPerFeatureQCMetrics}, further arguments to pass to \code{\link{perFeatureQCMetrics}}.
#'
#' @return
#' \code{x} is returned with the QC metrics added to the row or column metadata.
#' 
#' If possible the subsets are appended to the rowData of the \code{\link{SingleCellExperiment}} or \code{\link{SummarizedExperiment}}.
#'
#' @details
#' These functions are simply wrappers around \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}, respectively.
#' The computed QC metrics are automatically appended onto the existing \code{\link{colData}} or \code{\link{rowData}}.
#' No protection is provided against duplicated column names.
#'
#' \code{addPerCellQC} and \code{addPerFeatureQC} are exactly the same functions, \emph{sans} the \code{Metrics} at the end of their names.
#' They were added in the tempestuous youth of this package when naming was fast and loose.
#' These can be considered to be soft-deprecated in favor of the longer forms. 
#'
#' @author Aaron Lun, Lluís Revilla Sancho
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- addPerCellQCMetrics(example_sce,
#'                                    subsets = list(group1 = 1:5,
#'                                             group2 = c("Gene_0001", "Gene_2000"))
#' )
#' colData(example_sce)
#'
#' example_sce <- addPerFeatureQCMetrics(example_sce)
#' rowData(example_sce)
#' 
#' @seealso
#' \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}, which do the actual work.
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SummarizedExperiment rowData rowData<-
addPerCellQCMetrics <- function(x, subsets = NULL, ...) {
    colData(x) <- cbind(colData(x), perCellQCMetrics(x, subsets = subsets, ...))
    features <- featureSelected(x, subsets = subsets)
    if (ncol(features) > 0L) {
        row_data <- cbind(rowData(x), features)
        rownames(row_data) <- rownames(rowData(x))
        rowData(x) <- row_data
    }
    x
}

#' @importFrom S4Vectors make_zero_col_DFrame
featureSelected <- function(x, subsets) {
    n_features <- nrow(x)
    if (missing(subsets) || is.null(subsets)) {
        return(make_zero_col_DFrame(n_features))
    }
    
    features <- rownames(rowData(x))
    subsets_logical <- lapply(subsets, FUN = function(x, target) {
        if (length(x) == length(target) && is.logical(x) && !anyNA(x)) {
            return(x)
        }
        if (is.character(x)) {
            x <- which(x == features)
        }
        target[x] <- TRUE
        target
    }, target = vector("logical", n_features))
    DataFrame(subsets_logical)
}


#' @export
#' @rdname addPerCellQCMetrics
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment rowData rowData<-
addPerFeatureQCMetrics <- function(x, ...) {
    rowData(x) <- cbind(rowData(x), perFeatureQCMetrics(x, ...))
    x
}

#' @export
#' @rdname addPerCellQCMetrics
addPerCellQC <- addPerCellQCMetrics

#' @export
#' @rdname addPerCellQCMetrics
addPerFeatureQC <- addPerFeatureQCMetrics
