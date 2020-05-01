#' Sum counts across feature sets
#' 
#' Sum together expression values (by default, counts) for each feature set in each cell.
#'
#' @param x A numeric matrix of counts containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a count matrix.
#' @param ids A factor of length \code{nrow(x)}, specifying the set to which each feature in \code{x} belongs.
#'
#' Alternatively, a list of integer or character vectors, where each vector specifies the indices or names of features in a set.
#' @param average Logical scalar indicating whether the average should be computed instead of the sum.
#' @param subset.row An integer, logical or character vector specifying the features to use.
#' Defaults to all features.
#' @param subset.col An integer, logical or character vector specifying the cells to use.
#' Defaults to all cells with non-\code{NA} entries of \code{ids}.
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the matrix of counts
#' (or any other expression quantity that can be meaningfully summed).
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether summation should be parallelized.
#' @param ... For the \code{sumCountsAcrossFeatures} generic, further arguments to be passed to specific methods.
#' 
#' For the SummarizedExperiment method, further arguments to be passed to the ANY method.
#' @param subset_row,subset_col,exprs_values Soft-deprecated equivalents of the arguments described above.
#'
#' @return 
#' A count matrix is returned with one row per level of \code{ids}.
#' In each cell, counts for all features in the same set are summed together (or averaged, if \code{average=TRUE}).
#' Rows are ordered according to \code{levels(ids)}.
#'
#' @details
#' This function provides a convenient method for aggregating counts across multiple rows for each cell.
#' Several possible applications are listed below:
#' \itemize{
#' \item Using a list of genes in \code{ids}, we can obtain a summary expression value for all genes in one or more gene sets.
#' This allows the activity of various pathways to be compared across cells.
#' \item Genes with multiple mapping locations in the reference will often manifest as multiple rows with distinct Ensembl/Entrez IDs.
#' These counts can be aggregated into a single feature by setting the shared identifier (usually the gene symbol) as \code{ids}.
#' \item It is theoretically possible to aggregate transcript-level counts to gene-level counts with this function.
#' However, it is often better to do so with dedicated functions (e.g., from the \pkg{tximport} or \pkg{tximeta} packages) that account for differences in length across isoforms.
#' }
#'
#' The behaviour of this function is equivalent to that of \code{\link{rowsum}}.
#' However, this function can operate on any matrix representation in \code{object},
#' and can do so in a parallelized manner for large matrices without resorting to block processing.
#'
#' If \code{ids} is a factor, any \code{NA} values are implicitly ignored and will not be considered or reported.
#' This may be useful, e.g., to remove undesirable feature sets by setting their entries in \code{ids} to \code{NA}.
#'
#' Setting \code{average=TRUE} will compute the average in each set rather than the sum.
#' This is particularly useful if \code{x} contains expression values that have already been normalized in some manner,
#' as computing the average avoids another round of normalization to account for differences in the size of each set.
#'
#' @author Aaron Lun
#' @name sumCountsAcrossFeatures
#'
#' @seealso 
#' \code{\link{aggregateAcrossFeatures}}, to perform additional aggregation of row-level metadata.
#'
#' \code{\link{numDetectedAcrossFeatures}}, to compute the number of detected features per cell.
#'
#' @examples
#' example_sce <- mockSCE()
#' ids <- sample(LETTERS, nrow(example_sce), replace=TRUE)
#' out <- sumCountsAcrossFeatures(example_sce, ids)
#' str(out)
NULL

#' @importFrom BiocParallel SerialParam 
.sum_counts_across_features <- function(x, ids, subset.row=NULL, subset.col=NULL, 
    average=FALSE, BPPARAM=SerialParam(), subset_row=NULL, subset_col=NULL) 
{
    subset.row <- .replace(subset.row, subset_row)
    subset.col <- .replace(subset.col, subset_col)
    .sum_across_features(x, ids, subset.row=subset.row, subset.col=subset.col, average=average, BPPARAM=BPPARAM)
}

#' @importFrom BiocParallel SerialParam bplapply
.sum_across_features <- function(x, ids, subset.row=NULL, subset.col=NULL, 
    average=FALSE, BPPARAM=SerialParam(), modifier=NULL) 
{
    if (is.list(ids)) {
        ids <- lapply(ids, .subset2index, target=x, byrow=TRUE)
        runs <- lengths(ids)
        genes <- unlist(ids)
        names <- names(ids)
    } else {
        if (length(ids)!=nrow(x)) {
            stop("'ids' should be of length equal to 'nrow(x)'")            
        }
        ids <- factor(ids)
        genes <- order(ids, na.last=NA)
        runs <- as.integer(table(ids))
        names <- levels(ids)
    }

    if (!is.null(subset.row)) {
        subset.row <- .subset2index(subset.row, x, byrow=TRUE)
        keep <- genes %in% subset.row
        genes <- genes[keep]
        kept <- findInterval(which(keep), cumsum(runs), left.open=TRUE)
        runs <- tabulate(kept+1L, nbins=length(runs))
    }

    by.core <- .splitColsByWorkers(x, BPPARAM=BPPARAM, subset.col=subset.col)
    if (!is.null(modifier)) {
        by.core <- lapply(by.core, modifier)
    }
    out_list <- bplapply(by.core, FUN=sum_row_counts, genes=genes, runs=runs, BPPARAM=BPPARAM)

    out <- do.call(cbind, out_list)
    rownames(out) <- names
    colnames(out) <- unlist(lapply(by.core, colnames))
    if (average) {
        out <- out/runs
    }

    out
}

#' @export
#' @rdname sumCountsAcrossFeatures
setGeneric("sumCountsAcrossFeatures", function(x, ...) standardGeneric("sumCountsAcrossFeatures"))

#' @export
#' @rdname sumCountsAcrossFeatures
setMethod("sumCountsAcrossFeatures", "ANY", .sum_counts_across_features)

#' @export
#' @rdname sumCountsAcrossFeatures
#' @importFrom SummarizedExperiment assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod("sumCountsAcrossFeatures", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .sum_counts_across_features(assay(x, assay.type), ...)    
})
