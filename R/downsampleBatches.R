#' Downsample batches to equal coverage
#'
#' A convenience function to downsample all batches so that the average per-cell total count is the same across batches.
#' This mimics the downsampling functionality of \code{cellranger aggr}.
#' 
#' @param ... Two or more count matrices, where each matrix represents data from a separate batch.
#' Alternatively, a single filtered count matix containing cells from all batches, in which case \code{batch} should be specified.
#' @param batch A factor of length equal to the number of columns in the sole entry of \code{...},
#' specifying the batch of origin for each column of the matrix.
#' Ignored if there are multiple entries in \code{...}.
#' @param method String indicating how the average total should be computed.
#' The geometric mean is computed with a pseudo-count of 1.
#' @param bycol A logical scalar indicating whether downsampling should be performed on a column-by-column basis, 
#' see \code{?\link{downsampleMatrix}} for more details.
#'
#' @return
#' If \code{...} contains two or more matrices, a \linkS4class{List} of downsampled matrices is returned.
#'
#' Otherwise, if \code{...} contains only one matrix, the downsampled matrix is returned directly.
#'
#' @details
#' Downsampling batches with strong differences in sequencing coverage can make it easier to compare them to each other,
#' reducing the burden on the normalization and batch correction steps.
#' This is especially true when the number of cells cannot be easily controlled across batches,
#' resulting in large differences in per-cell coverage even when the total sequencing depth is the same.
#'
#' Generally speaking, the matrices in \code{...} should be filtered so that only libraries with cells are present.
#' This is most relevant to high-throughput scRNA-seq experiments (e.g., using  droplet-based protocols)
#' where the majority of libaries do not actually contain cells.
#' If these are not filtered out, downsampling will equalize coverage among the majority of empty libraries
#' rather than among cell-containing libraries.
#' 
#' @seealso
#' \code{\link{downsampleMatrix}}, which is called by this function under the hood.
#' 
#' @author Aaron Lun
#' 
#' @examples
#' sce1 <- mockSCE()
#' sce2 <- mockSCE()
#'
#' # Downsampling for multiple batches in a single matrix:
#' combined <- cbind(sce1, sce2)
#' batches <- rep(1:2, c(ncol(sce1), ncol(sce2)))
#' downsampled <- downsampleBatches(counts(combined), batch=batches)
#' downsampled[1:10,1:10]
#' 
#' # Downsampling for multiple matrices:
#' downsampled2 <- downsampleBatches(counts(sce1), counts(sce2))
#' downsampled2
#'
#' @export
#' @importFrom Matrix colSums
#' @importFrom S4Vectors List
#' @importFrom stats median
downsampleBatches <- function(..., batch=NULL, method=c("median", "mean", "geomean"), bycol=TRUE) {
    mats <- List(...)
    FUN <- switch(match.arg(method),
        mean=mean,
        median=median,
        geomean=function(x) { exp(mean(log(x+1))) }
    )

    if (length(mats)!=1L) {
        lib.sizes <- lapply(mats, colSums)
        measure <- vapply(lib.sizes, FUN, 0)
        down.prop <- min(measure)/measure

        for (i in seq_along(mats)) {
            mats[[i]] <- downsampleMatrix(mats[[i]], prop=down.prop[i], bycol=bycol)
        }
        mats

    } else {
        if (is.null(batch)) {
            stop("'batch' must be specified if '...' only contains one element")
        }
        lib.sizes <- colSums(mats[[1]]) 
        measure <- vapply(split(lib.sizes, batch), FUN, 0)
        down.prop <- (min(measure)/measure)[as.character(batch)]
        downsampleMatrix(mats[[1]], prop=down.prop)
    }
}
