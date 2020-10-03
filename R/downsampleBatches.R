#' Downsample batches to equal coverage
#'
#' A convenience function to downsample all batches so that the average per-cell total count is the same across batches.
#' This mimics the downsampling functionality of \code{cellranger aggr}.
#' 
#' @param ... Two or more count matrices, where each matrix has the same set of genes (rows)
#' and contains cells (columns) from a separate batch.
#'
#' Alternatively, one or more entries may be a \linkS4class{SummarizedExperiment}, 
#' in which case the count matrix is extracted from the assays according to \code{assay.type}.
#'
#' A list containing two or more of these matrices or SummarizedExperiments can also be supplied.
#'
#' Alternatively, a single count matrix or SummarizedExperiment can be supplied.
#' This is assumed to contain cells from all batches, in which case \code{batch} should also be specified.
#' @param batch A factor of length equal to the number of columns in the sole entry of \code{...},
#' specifying the batch of origin for each column of the matrix.
#' Ignored if there are multiple entries in \code{...}.
#' @param block If \code{...} contains multiple matrices or SummarizedExperiments,
#' this should be a character vector of length equal to the number of objects in \code{...},
#' specifying the blocking level for each object (see Details).
#'
#' Alternatively, if \code{...} contains a single object, 
#' this should be a character vector or factor of length specifing the blocking level for each cell in that object.
#' @param method String indicating how the average total should be computed.
#' The geometric mean is computed with a pseudo-count of 1.
#' @param bycol A logical scalar indicating whether downsampling should be performed on a column-by-column basis, 
#' see \code{?\link{downsampleMatrix}} for more details.
#' @param assay.type String or integer scalar specifying the assay of the SummarizedExperiment containing the count matrix,
#' if any SummarizedExperiments are present in \code{...}.
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
#' In more complex experiments,
#' batches can be organized into blocks where downsampling is performed to the lowest-coverage batch within each block.
#' This is most useful for larger datasets involving technical replicates for the same biological sample.
#' By setting \code{block=} to the biological sample, we can equalize coverage across replicates within each sample without forcing all samples to have the same coverage (e.g., to avoid loss of information if they are to be analyzed separately anyway).
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
#' @importFrom S4Vectors List
#' @importFrom stats median
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
downsampleBatches <- function(..., batch=NULL, block=NULL, method=c("median", "mean", "geomean"), bycol=TRUE, assay.type=1)
{
    mats <- .unpackLists(...)
    mats <- lapply(mats, FUN=function(x) {
        if (is(x, "SummarizedExperiment")) {
            x <- assay(x, assay.type)
        }
        x
    })

    FUN <- switch(match.arg(method),
        mean=mean,
        median=median,
        geomean=function(x) { exp(mean(log1p(x))) }
    )

    if (length(mats)!=1L) {
        if (is.null(block)) {
            List(.multi_object_downsample(mats, FUN=FUN, bycol=bycol))
        } else {
            if (length(block)!=length(mats)) {
                stop("'block' should be the same length as 'mats'")
            }
            by.block <- split(seq_along(mats), block)
            collected <- mats
            for (b in by.block) {
                collected[b] <- .multi_object_downsample(mats[b], FUN=FUN, bycol=bycol)
            }
            List(collected)
        }
    } else {
        if (is.null(batch)) {
            stop("'batch' must be specified if '...' only contains one element")
        }
        curmat <- mats[[1]]
        if (ncol(curmat)!=length(batch)) {
            stop("'length(batch)' must be equal to the number of cells") 
        }

        if (is.null(block)) {
            .single_object_downsample(curmat, batch, FUN=FUN, bycol=bycol)
        } else {
            if (length(block)!=ncol(curmat)) {
                stop("'length(block)' must be equal to the number of cells")
            }
            by.block <- split(seq_len(ncol(curmat)), block)
            collected <- by.block
            for (i in seq_along(by.block)) {
                b <- by.block[[i]]
                collected[[i]] <- .single_object_downsample(curmat[,b,drop=FALSE], batch[b], FUN=FUN, bycol=bycol)
            }
            output <- do.call(cbind, collected)
            output[,order(unlist(by.block)),drop=FALSE]    
        }
    }
}

#' @importFrom Matrix colSums
.multi_object_downsample <- function(mats, FUN, bycol) {
    lib.sizes <- lapply(mats, colSums)
    measure <- vapply(lib.sizes, FUN, 0)
    down.prop <- min(measure)/measure
    for (i in seq_along(mats)) {
        mats[[i]] <- downsampleMatrix(mats[[i]], prop=down.prop[i], bycol=bycol)
    }
    mats
}

.single_object_downsample <- function(mat, batch, FUN, bycol) { 
    by.batch <- split(seq_along(batch), batch)
    mats <- lapply(by.batch, function(i) mat[,i,drop=FALSE])
    mats <- .multi_object_downsample(mats, FUN=FUN, bycol=bycol)
    output <- do.call(cbind, mats)
    output[,order(unlist(by.batch)),drop=FALSE]
}
