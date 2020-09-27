#' Sum expression across groups of cells
#' 
#' Sum counts or average expression values for each feature across groups of cells.
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
#' Note that, prior to version 1.16.0, \code{sumCountsAcrossCells} would return a raw matrix.
#' This has now been wrapped in a \linkS4class{SummarizedExperiment} for consistency and to include per-group statistics.
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
#' attr(out, "ncells")
#'
#' batches <- sample(1:3, ncol(example_sce), replace=TRUE)
#' out2 <- sumCountsAcrossCells(example_sce, 
#'       DataFrame(label=ids, batch=batches))
#' head(out2)
#' attr(out2, "ids")
NULL

#' @importFrom BiocParallel SerialParam 
.sum_counts_across_cells <- function(x, ids, subset.row=NULL, subset.col=NULL,
    store.number="ncells", average=FALSE, BPPARAM=SerialParam(), 
    subset_row=NULL, subset_col=NULL, store_number=NULL)
{
    subset.row <- .replace(subset.row, subset_row)
    subset.col <- .replace(subset.col, subset_col)
    store.number <- .replace(store.number, store_number)

    .sum_across_cells_to_se(x, ids=ids, subset.row=subset.row, subset.col=subset.col,
        store.number=store.number, average=average, BPPARAM=BPPARAM)
}


#' @importFrom BiocParallel SerialParam 
#' @importFrom SummarizedExperiment SummarizedExperiment
.sum_across_cells_to_se <- function(x, ids, subset.row=NULL, subset.col=NULL,
    store.number="ncells", average=FALSE, BPPARAM=SerialParam(), modifier=NULL) 
{
    new.ids <- .process_ids(x, ids, subset.col)
    sum.out <- .sum_across_cells(x, ids=new.ids, subset.row=subset.row,
        average=average, BPPARAM=BPPARAM, modifier=modifier)

    mat <- sum.out$mat 
    mapping <- match(colnames(mat), as.character(new.ids))
    coldata <- .create_coldata(original.ids=ids, mapping=mapping, 
        freq=sum.out$freq, store.number=store.number)

    # Sync'ing the names as coldata is the source of truth.
    colnames(mat) <- rownames(coldata)

    output <- list(mat)
    names(output) <- sum.out$name
    SummarizedExperiment(output, colData=coldata)
}

.standardize_average_arg <- function(average) {
    if (isTRUE(average)) {
        average <- "mean"
    } else if (isFALSE(average)) {
        average <- "none"
    }
    match.arg(average, c("none", "mean", "median"))
}

#' @importFrom Matrix t 
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom BiocParallel SerialParam bplapply 
#' @importClassesFrom BiocParallel BiocParallelParam
.sum_across_cells <- function(x, ids, subset.row=NULL, average=FALSE, BPPARAM=SerialParam(), modifier=NULL) {
    average <- .standardize_average_arg(average)
    if (average=="median") {
        FUN <- .colmed
    } else {
        FUN <- .colsum
    }

    lost <- is.na(ids)
    subset.col <- if (any(lost)) which(!lost)
    by.core <- .splitRowsByWorkers(x, BPPARAM=BPPARAM, 
        subset.row=subset.row, subset.col=subset.col)

    if (!is.null(modifier)) { # used by numDetectedAcrossCells.
        by.core <- lapply(by.core, modifier)
    }

    # Avoid additional parallelization from DA methods.
    # We do the check primarily for testing purposes, as
    # flush tests with a FailParam defined in the GlobalEnv
    # don't get the FailParam class definition when this 
    # function is called via a SnowParam.
    oldBP <- getAutoBPPARAM()
    setAutoBPPARAM(SerialParam()) 
    if (is(oldBP, "BiocParallelParam")) {
        on.exit(setAutoBPPARAM(oldBP))
    }

    sub.ids <- ids[!lost]
    out <- bplapply(by.core, FUN=FUN, group=sub.ids, BPPARAM=BPPARAM)
    out <- do.call(rbind, out)

    freq <- table(sub.ids)
    freq <- as.integer(freq[colnames(out)])
    if (average=="mean") {
        out <- t(t(out)/freq)
    }

    if (average=="none") {
        name <- "sum" 
    } else { 
        name <- "average" 
    }

    list(mat=out, freq=freq, name=name)
}

##########################
##########################

#' @importFrom S4Vectors selfmatch
.df_to_factor <- function(ids) {
    o <- order(ids)
    x <- selfmatch(ids[o,,drop=FALSE]) 
    x[o] <- x
    x[Reduce("|", lapply(ids, is.na))] <- NA_integer_
    x
}

#' @importFrom methods is
.has_multi_ids <- function(ids) is(ids, "DataFrame")

.process_ids <- function(x, ids, subset.col) {    
    if (.has_multi_ids(ids)) {
        ids <- .df_to_factor(ids)
    } 
    if (ncol(x)!=length(ids)) {
        stop("length of 'ids' and 'ncol(x)' are not equal")
    }
    if (!is.null(subset.col)) {
        ids[!seq_along(ids) %in% .subset2index(subset.col, x, byrow=FALSE)] <- NA_integer_
    }
    ids
}

#' @importFrom S4Vectors DataFrame
.create_coldata <- function(original.ids, mapping, freq, store.number) {
    if (.has_multi_ids(original.ids)) {
        coldata <- original.ids[mapping,,drop=FALSE]
        rownames(coldata) <- NULL
    } else {
        coldata <- DataFrame(ids=original.ids[mapping])
        rownames(coldata) <- coldata$ids
    }

    if (!is.null(store.number)) {
        coldata[[store.number]] <- freq
    }
    coldata
}

##########################
##########################

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

##########################
##########################

setGeneric(".colsum", function(x, group) standardGeneric(".colsum"))

#' @importFrom Matrix rowSums
setMethod(".colsum", "ANY", function(x, group) {
    by.group <- split(seq_len(ncol(x)), group)
    out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE]))
    do.call(cbind, out)
})

#' @importFrom DelayedArray colsum
setMethod(".colsum", "matrix", function(x, group) {
    colsum(x, group)
})

#' @importFrom DelayedArray colsum
setMethod(".colsum", "DelayedMatrix", function(x, group) {
    colsum(x, group)
})

#' @importFrom DelayedMatrixStats rowMedians
.colmed <- function(x, group) {
    by.group <- split(seq_along(group), group)
    output <- matrix(0, nrow(x), length(by.group),
        dimnames=list(rownames(x), names(by.group)))

    for (i in names(by.group)) {
        current <- x[,by.group[[i]],drop=FALSE]
        current <- DelayedArray(current)
        output[,i] <- rowMedians(current)
    }
    output
}
