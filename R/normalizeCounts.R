#' Compute normalized expression values
#'
#' Compute (log-)normalized expression values by dividing counts for each cell by the corresponding size factor.
#'
#' @param x A numeric matrix-like object containing counts for cells in the columns and features in the rows.
#'
#' Alternatively, a \linkS4class{SingleCellExperiment} or \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the count matrix.
#' @param size.factors A numeric vector of cell-specific size factors.
#' Alternatively \code{NULL}, in which case the size factors are extracted or computed from \code{x}.
#' @param log Logical scalar indicating whether normalized values should be log2-transformed.
#' @param pseudo.count Numeric scalar specifying the pseudo-count to add when log-transforming expression values.
#' @param center.size.factors Logical scalar indicating whether size factors should be centered at unity before being used.
#' @param subset.row A vector specifying the subset of rows of \code{x} for which to return a result.
#' @param downsample Logical scalar indicating whether downsampling should be performed prior to scaling and log-transformation.
#' @param down.target Numeric scalar specifying the downsampling target when \code{downsample=TRUE}.
#' If \code{NULL}, this is defined by \code{down.prop} and a warning is emitted.
#' @param down.prop Numeric scalar between 0 and 1 indicating the quantile to use to define the downsampling target.
#' Only used when \code{downsample=TRUE}.
#' @param ... For the generic, arguments to pass to specific methods.
#'
#' For the SummarizedExperiment method, further arguments to pass to the ANY or \linkS4class{DelayedMatrix} methods.
#' 
#' For the SingleCellExperiment method, further arguments to pass to the SummarizedExperiment method.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how library size factor calculations should be parallelized.
#' Only used if \code{size.factors} is not specified.
#' @param exprs_values,size_factors,pseudo_count,center_size_factors,subset_row,down_target,down_prop
#' Soft-deprecated equivalents to the arguments described previously.
#'
#' @details 
#' Normalized expression values are computed by dividing the counts for each cell by the size factor for that cell.
#' This removes cell-specific scaling biases due to differences in sequencing coverage, capture efficiency or total RNA content.
#' If \code{log=TRUE}, log-normalized values are calculated by adding \code{pseudo.count} to the normalized count and performing a log2-transformation.
#'
#' If no size factors are supplied, they are determined automatically from \code{x}:
#' \itemize{
#' \item For count matrices and \linkS4class{SummarizedExperiment} inputs,
#' the sum of counts for each cell is used to compute a size factor via the \code{\link{librarySizeFactors}} function.
#' \item For \linkS4class{SingleCellExperiment} instances, the function searches for \code{\link{sizeFactors}} from \code{x}.
#' If none are available, it defaults to library size-derived size factors.
#' }
#' If \code{size.factors} are supplied, they will override any size factors present in \code{x}.
#'
#' @section Centering the size factors:
#' If \code{center.size.factors=TRUE}, size factors are centred at unity prior to calculation of normalized expression values.
#' This ensures that the computed expression values can be interpreted as being on the same scale as original counts.
#' We can then compare abundances between features normalized with different sets of size factors; the most common use of this fact is in the comparison between spike-in and endogenous abundances when modelling technical noise (see \code{\link[scran]{modelGeneVarWithSpikes}} package for an example).
#'
#' More generally, when \code{log=TRUE}, centering of the size factors ensures that the value of \code{pseudo.count} can be interpreted as being on the same scale as the counts, i.e., the pseudo-count can actually be thought of as a \emph{count}.
#' This is important as it implies that the pseudo-count's impact will diminish as sequencing coverage improves.
#' Thus, if the size factors are centered, differences between log-normalized expression values will more closely approximate the true log-fold change with increasing coverage, whereas this would not be true of other metrics like log-CPMs with a fixed offset.
#'
#' The disadvantage of using centered size factors is that the expression values are not directly comparable across different calls to \code{\link{normalizeCounts}}, typically for multiple batches.
#' In theory, this is not a problem for metrics like the CPM, but in practice, we have to apply batch correction methods anyway to perform any joint analysis - see \code{\link[batchelor]{multiBatchNorm}} for more details. 
#'
#' @section Downsampling instead of scaling:
#' If \code{downsample=TRUE}, counts for each cell are randomly downsampled instead of being scaled.
#' This is occasionally useful for avoiding artifacts caused by scaling count data with a strong mean-variance relationship.
#' Each cell is downsampled according to the ratio between \code{down.target} and that cell's size factor.
#' (Cells with size factors below the target are not downsampled and are directly scaled by this ratio.)
#' If \code{log=TRUE}, a log-transformation is also performed after adding \code{pseudo.count} to the downsampled counts.
#'
#' We automatically set \code{down.target} to the 1st percentile of size factors across all cells involved in the analysis,
#' but this is only appropriate if the resulting expression values are not compared across different \code{normalizeCounts} calls.
#' To obtain expression values that are comparable across different \code{normalizeCounts} calls
#' (e.g., in \code{\link[scran]{modelGeneVarWithSpikes}} or \code{\link[batchelor]{multiBatchNorm}}),
#' \code{down_target} should be manually set to a constant target value that can be considered a low size factor in every call.
#'  
#' @return A numeric matrix-like object containing (log-)normalized expression values.
#' This has the same dimensions as \code{x} (unless \code{subset.row} is specified)
#' and is of the same class as the original count matrix.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{logNormCounts}}, which wraps this function for convenient use with SingleCellExperiment instances.
#'
#' \code{\link{librarySizeFactors}}, to compute the default size factors.
#'
#' \code{\link[DropletUtils]{downsampleMatrix}}, to perform the downsampling.
#' @examples
#' example_sce <- mockSCE()
#'
#' # Standard scaling + log-transformation:
#' normed <- normalizeCounts(example_sce)
#' normed[1:5,1:5]
#'
#' # Scaling without transformation:
#' normed <- normalizeCounts(example_sce, log=FALSE)
#' normed[1:5,1:5]
#'
#' # Downscaling with transformation:
#' normed <- normalizeCounts(example_sce, downsample=TRUE)
#' normed[1:5,1:5]
#'
#' # Using custom size factors:
#' with.meds <- computeMedianFactors(example_sce)
#' normed <- normalizeCounts(with.meds)
#' normed[1:5,1:5]
#' 
#' @name normalizeCounts
NULL

#' @export
#' @rdname normalizeCounts
setGeneric("normalizeCounts", function(x, ...) standardGeneric("normalizeCounts"))

#' @export
#' @rdname normalizeCounts
#' @importFrom BiocParallel SerialParam
setMethod("normalizeCounts", "ANY", function(x, size.factors=NULL, 
    log=TRUE, pseudo.count=1, center.size.factors=TRUE, subset.row=NULL,
    downsample=FALSE, down.target=NULL, down.prop=0.01, BPPARAM=SerialParam(),
    size_factors=NULL, pseudo_count=NULL, center_size_factors=NULL,
    subset_row=NULL, down_target=NULL, down_prop=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    size.factors <- .replace(size.factors, size_factors)
    center.size.factors <- .replace(center.size.factors, center_size_factors)
    down.target <- .replace(down.target, down_target)
    down.prop <- .replace(down.prop, down_prop)
    pseudo.count <- .replace(pseudo.count, pseudo_count)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    if (nrow(x)==0L) {
        return(x + 0) # coerce to numeric.
    }

    size.factors <- .get_default_sizes(x, size.factors, center.size.factors, BPPARAM=BPPARAM)
    if (length(size.factors)!=ncol(x)) {
        stop("number of size factors does not equal 'ncol(x)'")
    }
    if (!all(is.finite(size.factors) & size.factors > 0)) {
        stop("size factors should be positive")
    }

    if (downsample) {
        down.out <- .downsample_counts(x, size.factors, down.prop=down.prop, down.target=down.target)
        x <- down.out$x
        size.factors <- down.out$size.factors
    }

    .internal_transformer(x, size.factors, log, pseudo.count) 
})

.get_default_sizes <- function(x, size.factors, center.size.factors, ...) {
    if (is.null(size.factors)) {
        size.factors <- librarySizeFactors(x, ...)
    }
    .center.size.factors(size.factors, center.size.factors)
}

.center.size.factors <- function(size.factors, center.size.factors) {
    if (center.size.factors) {
        size.factors <- size.factors/mean(size.factors)
    }
    size.factors
}

#' @importFrom stats quantile
.downsample_counts <- function(x, size.factors, down.prop, down.target) {
    if (is.null(down.target)) {
        down.target <- quantile(size.factors, probs=down.prop)
        warning("'down.target' defined as the 1st percentile of size factors")
    }
    down_rate <- pmin(1, down.target/size.factors)
    x <- DropletUtils::downsampleMatrix(x, down_rate, bycol=TRUE)
    size.factors <- size.factors * down_rate/down.target
    list(x=x, size.factors=size.factors)
}

###########################################

#' @export
#' @rdname normalizeCounts
#' @importFrom SummarizedExperiment assay 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("normalizeCounts", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    normalizeCounts(assay(x, assay.type), ...)
})

#' @export
#' @rdname normalizeCounts
#' @importFrom BiocGenerics sizeFactors
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod("normalizeCounts", "SingleCellExperiment", function(x, size.factors=NULL, ...) {
    if (is.null(size.factors)) {
        size.factors <- sizeFactors(x)
    }
    callNextMethod(x=x, size.factors=size.factors, ...)
})

###########################################

setGeneric(".internal_transformer", function(x, ...) standardGeneric(".internal_transformer"))

#' @importFrom Matrix t
setMethod(".internal_transformer", "ANY", function(x, size.factors, log, pseudo.count) {
    norm_exprs <- t(t(x) / size.factors)
    if (log) {
        norm_exprs <- log2(norm_exprs + pseudo.count)
    }
    norm_exprs
})

#' @importClassesFrom Matrix dgTMatrix
setMethod(".internal_transformer", "dgCMatrix", function(x, size.factors, log, pseudo.count) {
    if (log && pseudo.count!=1) {
        callNextMethod()
    } else {
        .transform_sparse(x, rep(size.factors, diff(x@p)), log, pseudo.count)
    }
})

#' @importClassesFrom Matrix dgTMatrix
setMethod(".internal_transformer", "dgTMatrix", function(x, size.factors, log, pseudo.count) {
    if (log && pseudo.count!=1) {
        callNextMethod()
    } else {
        .transform_sparse(x, size.factors[x@j+1L], log, pseudo.count)
    }
})

.transform_sparse <- function(x, expanded_sf, log, pseudo.count) {
    x@x <- x@x/expanded_sf
    if (log) {
        x@x <- log2(x@x + pseudo.count)
    }
    x
}
