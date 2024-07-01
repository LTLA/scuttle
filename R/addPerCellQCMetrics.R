#' Add QC metrics to a SummarizedExperiment
#'
#' Convenient utilities to compute QC metrics and add them to a \linkS4class{SummarizedExperiment}'s row or column metadata.
#'
#' @param x A \linkS4class{SummarizedExperiment} object or one of its subclasses.
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
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' example_sce <- addPerCellQCMetrics(example_sce)
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
addPerCellQCMetrics <- function(x, ...) {
    colData(x) <- cbind(colData(x), perCellQCMetrics(x, ...))
    features <- featureSelected(x, ...)
    if (!is.null(features)) {
      rowData(x) <- cbind(rowData(x), features)
    }
    x
}

#' @export
#' @rdname featureSelected
setGeneric("featureSelected", function(x, ...) standardGeneric("featureSelected"))

#' @export
#' @rdname featureSelected
setMethod("featureSelected", "SummarizedExperiment", function(x, ..., subsets) {
  if (missing(subsets) || is.null(subsets)) {
    return(NULL)
  }
  dummy <- vector("logical", nrow(x))
  
  subsets_logical <- lapply(subsets, FUN = function(x, target) {
    target[x] <- TRUE
    target
  }, target = dummy)
  DataFrame(simplify2array(subsets_logical))
})


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
