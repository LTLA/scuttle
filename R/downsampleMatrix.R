#' Downsample a count matrix
#'
#' Downsample a count matrix to a desired proportion, either on a whole-matrix or per-cell basis.
#' 
#' @param x An integer or numeric matrix-like object containing counts.
#' @param prop A numeric scalar or, if \code{bycol=TRUE}, a vector of length \code{ncol(x)}.
#' All values should lie in [0, 1] specifying the downsampling proportion for the matrix or for each cell.
#' @param bycol A logical scalar indicating whether downsampling should be performed on a column-by-column basis.
#' @param sink A \linkS4class{RealizationSink} object specifying the format of the downsampled matrix should be returned.
#' 
#' @return An numeric matrix-like object of downsampled counts.
#' This is usually of the same class as \code{x}, unless \code{sink} is set in which case it is a \linkS4class{DelayedMatrix}.
#' 
#' @details
#' Given multiple batches of very different sequencing depths, 
#' it can be beneficial to downsample the deepest batches to match the coverage of the shallowest batches. 
#' This avoids differences in technical noise that can drive clustering by batch.
#' 
#' If \code{bycol=TRUE}, sampling without replacement is performed on the count vector for each cell.
#' This yields a new count vector where the total is equal to \code{prop} times the original total count. 
#' Each count in the returned matrix is guaranteed to be smaller than the original value in \code{x}.
#' Different proportions can be specified for different cells by setting \code{prop} to a vector;
#' in this manner, downsampling can be used as an alternative to scaling for per-cell normalization.
#' 
#' If \code{bycol=FALSE}, downsampling without replacement is performed on the entire matrix.
#' This yields a new matrix where the total count across all cells is equal to \code{prop} times the original total.
#' The new total count for each cell may not be exactly equal to \code{prop} times the original value,
#' which may or may not be more appropriate than \code{bycol=TRUE} for particular applications.
#' 
#' By default, the output format will be a dense matrix for an ordinary matrix \code{x},
#' and a sparse matrix for any sparse format (i.e., if \code{is_sparse(x)} returns \code{TRUE}).
#' This can be overridden to specify custom formats with \code{sink}, e.g., for HDF5-backed matrices,
#' which may be helpful for dealing with large matrices where the downsampled result may not fit in memory.
#' 
#' Note that this function was originally implemented in the \pkg{scater} package as \code{downsampleCounts},
#' was moved to the \pkg{DropletUtils} package as \code{downsampleMatrix},
#' and finally found a home here.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link[DropletUtils]{downsampleReads}} in the \pkg{DropletUtils} package, 
#' which downsamples reads rather than observed counts.
#'
#' \code{\link{normalizeCounts}}, where downsampling can be used as an alternative to scaling normalization.
#' 
#' @examples
#' sce <- mockSCE()
#' sum(counts(sce))
#' 
#' downsampled <- downsampleMatrix(counts(sce), prop = 0.5, bycol=FALSE)
#' sum(downsampled)
#'
#' downsampled2 <- downsampleMatrix(counts(sce), prop = 0.5, bycol=TRUE)
#' sum(downsampled2)
#'
#' @export
#' @importFrom DelayedArray makeNindexFromArrayViewport write_block currentViewport
#' @importClassesFrom DelayedArray DelayedArray
downsampleMatrix <- function(x, prop, bycol=TRUE, sink=NULL) {
    if (bycol) {
        prop <- rep(prop, length.out = ncol(x))
        FUN <- function(block) {
            vp <- currentViewport()
            cols <- makeNindexFromArrayViewport(vp, expand.RangeNSBS=TRUE)[[2]]
            if (!is.null(cols)) {
                prop <- prop[cols]
            }

            downed <- downsample_column(block, prop)

            if (!is.null(sink)) {
                write_block(sink, vp, as.matrix(downed))
            } else {
                downed
            }
        }
    } else {
        # Mimic the downsample_vector_part constructor 
        # in inst/include/scuttle/downsample_vector.h 
        total <- round(sum(x))
        required <- round(total * min(1, prop))

        FUN <- function(block) {
            out <- downsample_matrix(block, total, required)
            downed <- out[[1]]
            total <<- total - out[[2]]
            required <<- required - sum(out[[1]])

            if (!is.null(sink)) {
                vp <- currentViewport()
                write_block(sink, vp, as.matrix(downed))
                NULL
            } else {
                downed
            }
        }
    }

    # This MUST be excuted in serial, given the <<- and the correct behavior of
    # the random stream (can't be bothered to use dqrng here).
    out <- colBlockApply(x, FUN=FUN, BPPARAM=NULL) 

    if (is.null(sink)) {
        if (length(out)!=1L) {
            out <- do.call(cbind, out)
        } else {
            out <- out[[1]]
        }
    } else {
        out <- as(sink, "DelayedArray")
    }

    if (!is.null(out)) {
        dimnames(out) <- dimnames(x)
    }
    out
}
