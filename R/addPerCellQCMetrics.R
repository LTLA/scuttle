#' Add QC metrics to a SummarizedExperiment
#'
#' Convenient utilities to compute QC metrics and add them to a \linkS4class{SummarizedExperiment}'s row or column metadata.
#'
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
#' @inheritParams perCellQCMetrics
#' @param subset.prefix String containing the prefix for the names of the columns of \code{\link{rowData}} that specify which genes belong to each subset.
#' If \code{NULL}, these subset identity columns are not added to the \code{\link{rowData}}.
#' @param ... For \code{addPerCellQCMetrics}, further arguments to pass to \code{\link{perCellQCMetrics}}.
#' 
#' For \code{addPerFeatureQCMetrics}, further arguments to pass to \code{\link{perFeatureQCMetrics}}.
#'
#' @return
#' \code{x} is returned with the QC metrics added to the row or column metadata.
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
#' example_sce <- addPerCellQCMetrics(
#'      example_sce,
#'      subsets = list(group1 = 1:5, group2 = c("Gene_0001", "Gene_2000"))
#' )
#' colData(example_sce)
#' rowData(example_sce)
#'
#' example_sce <- addPerFeatureQCMetrics(example_sce)
#' rowData(example_sce)
#' 
#' @seealso
#' \code{\link{perCellQCMetrics}} and \code{\link{perFeatureQCMetrics}}, which do the actual work.
#'
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment colData colData<- rowData rowData<-
addPerCellQCMetrics <- function(x, subsets = NULL, ..., subset.prefix = "subset_") {
    colData(x) <- cbind(colData(x), perCellQCMetrics(x, subsets = subsets, ...))

    if (!is.null(subset.prefix) && length(subsets)) {
        row_data <- rowData(x)
        for (s in names(subsets)) {
            in.subset <- logical(nrow(x))
            in.subset[.subset2index(subsets[[s]], x)] <- TRUE
            row_data[[paste0(subset.prefix, s)]] <- in.subset
        }
        rowData(x) <- row_data
    }

    x
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
