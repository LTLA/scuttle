#############################
# Metric calculator generics.

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
