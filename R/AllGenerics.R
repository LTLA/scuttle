#############################
# Metric calculator generics.

#' @export
#' @rdname perCellQCMetrics
setGeneric("perCellQCMetrics", function(x, ...) standardGeneric("perCellQCMetrics"))

#' @export
#' @rdname perFeatureQCMetrics
setGeneric("perFeatureQCMetrics", function(x, ...) standardGeneric("perFeatureQCMetrics"))

#' @export
#' @rdname getVarianceExplained
setGeneric("getVarianceExplained", function(x, ...) standardGeneric("getVarianceExplained"))

#' @export
#' @rdname nexprs
setGeneric("nexprs", function(x, ...) standardGeneric("nexprs"))

####################################
# Dimensionality reduction generics.

#' @export
#' @rdname runPCA
setGeneric("calculatePCA", function(x, ...) standardGeneric("calculatePCA"))

#' @export
#' @rdname runTSNE
setGeneric("calculateTSNE", function(x, ...) standardGeneric("calculateTSNE"))

#' @export
#' @rdname runUMAP
setGeneric("calculateUMAP", function(x, ...) standardGeneric("calculateUMAP"))

#' @export
#' @rdname runMDS
setGeneric("calculateMDS", function(x, ...) standardGeneric("calculateMDS"))

#' @export
#' @rdname runNMF
setGeneric("calculateNMF", function(x, ...) standardGeneric("calculateNMF"))

#' @export
#' @rdname runDiffusionMap
setGeneric("calculateDiffusionMap", function(x, ...) standardGeneric("calculateDiffusionMap"))

#################################
# Aggregation-related generics.

#' @export
#' @rdname sumCountsAcrossFeatures
setGeneric("sumCountsAcrossFeatures", function(x, ...) standardGeneric("sumCountsAcrossFeatures"))

#' @export
#' @rdname sumCountsAcrossCells
setGeneric("sumCountsAcrossCells", function(x, ...) standardGeneric("sumCountsAcrossCells"))

#' @export
#' @rdname sumCountsAcrossCells
setGeneric("aggregateAcrossCells", function(x, ...) standardGeneric("aggregateAcrossCells"))

#' @export
#' @rdname numDetectedAcrossCells
setGeneric("numDetectedAcrossCells", function(x, ...) standardGeneric("numDetectedAcrossCells"))

#' @export
#' @rdname numDetectedAcrossFeatures
setGeneric("numDetectedAcrossFeatures", function(x, ...) standardGeneric("numDetectedAcrossFeatures"))
