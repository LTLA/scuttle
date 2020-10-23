#' Sum expression across groups of cells
#' 
#' Sum counts or average expression values for each feature across groups of cells.
#' This function is deprecated; use \code{\link{summarizeAssayByGroup}} instead.
#'
#' @param x A numeric matrix of expression values (usually counts) containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param ids A factor specifying the group to which each cell in \code{x} belongs.
#' Alternatively, a \linkS4class{DataFrame} of such vectors or factors, 
#' in which case each unique combination of levels defines a group. 
#' @param subset.row An integer, logical or character vector specifying the features to use.
#' If \code{NULL}, defaults to all features.
#' For the \linkS4class{SingleCellExperiment} method, this argument will not affect alternative Experiments,
#' where aggregation is always performed for all features (or not at all, depending on \code{use_alt_exps}).
#' @param subset.col An integer, logical or character vector specifying the cells to use.
#' Defaults to all cells with non-\code{NA} entries of \code{ids}.
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the matrix of counts
#' (or any other expression quantity that can be meaningfully summed).
#' @param average Logical scalar indicating whether the average should be computed instead of the sum.
#' Alternatively, a string containing \code{"mean"}, \code{"median"} or \code{"none"}, specifying the type of average.
#' (\code{"none"} is equivalent to \code{FALSE}.)
#' @param store.number String specifying the field of the output \code{\link{colData}} to store the number of cells in each group.
#' If \code{NULL}, nothing is stored.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether summation should be parallelized.
#' @param ... For the generics, further arguments to be passed to specific methods.
#' 
#' For the SummarizedExperiment method, further arguments to be passed to the ANY method.
#' @param subset_row,subset_col,exprs_values,store_number Soft-deprecated equivalents to the arguments described above.
#'
#' @return 
#' A SummarizedExperiment is returned with one column per level of \code{ids}.
#' Each entry of the assay contains the sum or average across all cells in a given group (column) for a given feature (row).
#' Columns are ordered by \code{levels(ids)} and the number of cells per level is reported in the \code{"ncells"} column metadata.
#' For DataFrame \code{ids}, each column corresponds to a unique combination of levels (recorded in the \code{\link{colData}}).
#'
#' @details
#' These functions provide a convenient method for summing or averaging expression values across multiple columns for each feature.
#' A typical application would be to sum counts across all cells in each cluster to obtain \dQuote{pseudo-bulk} samples for further analyses, e.g., differential expression analyses between conditions.
#'
#' The behaviour of \code{sumCountsAcrossCells} is equivalent to that of \code{\link{colsum}}.
#' However, this function can operate on any matrix representation in \code{object};
#' can do so in a parallelized manner for large matrices without resorting to block processing;
#' and can natively support combinations of multiple factors in \code{ids}.
#'
#' Any \code{NA} values in \code{ids} are implicitly ignored and will not be considered during summation.
#' This may be useful for removing undesirable cells by setting their entries in \code{ids} to \code{NA}.
#' Alternatively, we can explicitly select the cells of interest with \code{subset_col}.
#' 
#' Setting \code{average=TRUE} will compute the average in each set rather than the sum.
#' This is particularly useful if \code{x} contains expression values that have already been normalized in some manner,
#' as computing the average avoids another round of normalization to account for differences in the size of each set.
#' The same effect is obtained by setting \code{average="mean"},
#' while setting \code{average="median"} will instead compute the median across all cells.
#'
#' @author Aaron Lun
#' @name sumCountsAcrossCells
#'
#' @seealso
#' \code{\link{aggregateAcrossCells}}, which also combines information in the \code{colData}.
#'
#' \code{\link{numDetectedAcrossCells}}, which computes the number of cells with detected expression in each group.
#'
#' @examples
#' example_sce <- mockSCE()
#' ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)
#'
#' out <- sumCountsAcrossCells(example_sce, ids)
#' head(out)
#'
#' batches <- sample(1:3, ncol(example_sce), replace=TRUE)
#' out2 <- sumCountsAcrossCells(example_sce, 
#'       DataFrame(label=ids, batch=batches))
#' head(out2)
NULL

#' @importFrom BiocParallel SerialParam 
#' @importFrom SummarizedExperiment assayNames<-
.sum_counts_across_cells <- function(x, ids, subset.row=NULL, subset.col=NULL,
    store.number="ncells", average=FALSE, BPPARAM=SerialParam(), 
    subset_row=NULL, subset_col=NULL, store_number=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    subset.col <- .replace(subset.col, subset_col)
    store.number <- .replace(store.number, store_number)

    average <- .average2statistic(average)
    output <- summarizeAssayByGroup(x, ids, subset.row=subset.row, subset.col=subset.col,
        statistics=average, store.number=store.number, BPPARAM=BPPARAM)

    if (average!="sum") {
        assayNames(output) <- "average"
    }
    output
}

.average2statistic <- function(average) {
    if (isTRUE(average)) {
        "mean"
    } else if (isFALSE(average)) {
        "sum"
    } else {
        average <- match.arg(average, c("none", "mean", "median"))
        if (average=="none") average <- "sum"
        average
    }
}

#' @export
#' @rdname sumCountsAcrossCells
setGeneric("sumCountsAcrossCells", function(x, ...) standardGeneric("sumCountsAcrossCells"))

#' @export
#' @rdname sumCountsAcrossCells
#' @importFrom BiocParallel SerialParam
setMethod("sumCountsAcrossCells", "ANY", .sum_counts_across_cells)

#' @export
#' @rdname sumCountsAcrossCells
#' @importFrom SummarizedExperiment assay
setMethod("sumCountsAcrossCells", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .sum_counts_across_cells(assay(x, assay.type), ...)
})
