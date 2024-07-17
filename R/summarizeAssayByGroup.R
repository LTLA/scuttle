#' Summarize an assay by group 
#' 
#' From an assay matrix, compute summary statistics for groups of cells.
#' A typical example would be to compute various summary statistics for clusters.
#'
#' @param x A numeric matrix containing features in rows and cells in columns.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#' @param ids A factor (or vector coercible into a factor) specifying the group to which each cell in \code{x} belongs.
#' Alternatively, a \linkS4class{DataFrame} of such vectors or factors, 
#' in which case each unique combination of levels defines a group. 
#' @param subset.row An integer, logical or character vector specifying the features to use.
#' If \code{NULL}, defaults to all features.
#' @param subset.col An integer, logical or character vector specifying the cells to use.
#' Defaults to all cells with non-\code{NA} entries of \code{ids}.
#' @param assay.type A string or integer scalar specifying the assay of \code{x} containing the assay to be summarized.
#' @param store.number String specifying the field of the output \code{\link{colData}} to store the number of cells in each group.
#' If \code{NULL}, nothing is stored.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether summation should be parallelized.
#' @param ... For the generics, further arguments to be passed to specific methods.
#' 
#' For the SummarizedExperiment method, further arguments to be passed to the ANY method.
#' @param statistics Character vector specifying the type of statistics to be computed, see Details.
#' @param threshold A numeric scalar specifying the threshold above which a gene is considered to be detected.
#'
#' @return 
#' A SummarizedExperiment is returned with one column per level of \code{ids}.
#' Each entry of the assay contains the sum or average across all cells in a given group (column) for a given feature (row).
#' Columns are ordered by \code{levels(ids)} and the number of cells per level is reported in the \code{"ncells"} column metadata.
#' For DataFrame \code{ids}, each column corresponds to a unique combination of levels (recorded in the \code{\link{colData}}).
#'
#' @details
#' These functions provide a convenient method for summing or averaging expression values across multiple columns for each feature.
#' A typical application would be to sum counts across all cells in each cluster to obtain \dQuote{pseudo-bulk} samples for further analyses, 
#' e.g., differential expression analyses between conditions.
#'
#' For each feature, the chosen assay can be aggregated by:
#' \itemize{
#' \item \code{"sum"}, the sum of all values in each group.
#' This makes the most sense for raw counts, to allow models to account for the mean-variance relationship.
#' \item \code{"mean"}, the mean of all values in each group.
#' This makes the most sense for normalized and/or transformed assays.
#' \item \code{"median"}, the median of all values in each group.
#' This makes the most sense for normalized and/or transformed assays, 
#' usually generated from large counts where discreteness is less of an issue.
#' \item \code{"num.detected"} and \code{"prop.detected"}, the number and proportion of values in each group that are non-zero.
#' This makes the most sense for raw counts or sparsity-preserving transformations.
#' }
#'
#' Any \code{NA} values in \code{ids} are implicitly ignored and will not be considered during summation.
#' This may be useful for removing undesirable cells by setting their entries in \code{ids} to \code{NA}.
#' Alternatively, we can explicitly select the cells of interest with \code{subset_col}.
#'
#' If \code{ids} is a factor and contains unused levels, they will not be represented as columns in the output.
#' 
#' @author Aaron Lun
#' @name summarizeAssayByGroup
#'
#' @seealso
#' \code{\link{aggregateAcrossCells}}, which also combines information in the \code{\link{colData}} of \code{x}.
#'
#' @examples
#' example_sce <- mockSCE()
#' ids <- sample(LETTERS[1:5], ncol(example_sce), replace=TRUE)
#'
#' out <- summarizeAssayByGroup(example_sce, ids)
#' out
#'
#' batches <- sample(1:3, ncol(example_sce), replace=TRUE)
#' out2 <- summarizeAssayByGroup(example_sce, 
#'       DataFrame(label=ids, batch=batches))
#' head(out2)
NULL

##########################
##########################

#' @importFrom BiocParallel SerialParam 
#' @importFrom SummarizedExperiment SummarizedExperiment
.summarize_assay_by_group <- function(x, ids, subset.row=NULL, subset.col=NULL,
    statistics=c("mean", "sum", "num.detected", "prop.detected", "median"),
    store.number="ncells", threshold=0, BPPARAM=SerialParam()) 
{
    new.ids <- .process_ids(x, ids, subset.col)
    sum.out <- .summarize_assay(x, ids=new.ids, subset.row=subset.row,
        statistics=match.arg(statistics, several.ok=TRUE),
        threshold=threshold, BPPARAM=BPPARAM)

    mat.out <- sum.out$summary
    mapping <- match(colnames(mat.out[[1]]), as.character(new.ids))
    coldata <- .create_coldata(original.ids=ids, mapping=mapping, 
        freq=sum.out$freq, store.number=store.number)

    # Sync'ing the names as coldata is the source of truth.
    for (i in seq_along(mat.out)) {
        colnames(mat.out[[i]]) <- rownames(coldata)
    }

    SummarizedExperiment(mat.out, colData=coldata)
}

#' @importFrom BiocParallel SerialParam 
#' @importFrom beachmat rowBlockApply
.summarize_assay <- function(x, ids, statistics, threshold=0, subset.row=NULL, BPPARAM=SerialParam()) {
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    if (is(x, "dgCMatrix") && all(statistics %in% c("sum", "mean", "num.detected", "prop.detected"))) {
        # Fast path hack to get around DelayedArray's inefficiencies.
        collected <- list()
        f <- factor(ids) # automatically drops unused levels.
        g <- as.integer(f) - 1L
        ngroups <- nlevels(f)
        freq <- tabulate(g + 1L, ngroups)
        names(freq) <- levels(f)

        if ("sum" %in% statistics || "mean" %in% statistics) {
            summed <- sparse_aggregate_sum(x@x, x@i, x@p, g, ngroups, nrow(x))
            if ("sum" %in% statistics) {
                collected$sum <- summed
            }
            if ("mean" %in% statistics) {
                collected$mean <- t(t(summed) / freq)
            }
        }

        if ("num.detected" %in% statistics || "prop.detected" %in% statistics) {
            detected <- sparse_aggregate_detected(x@x, x@i, x@p, g, ngroups, nrow(x))
            if ("num.detected" %in% statistics) {
                collected$num.detected <- detected
            }
            if ("prop.detected" %in% statistics) {
                collected$prop.detected <- t(t(detected) / freq)
            }
        }

        for (i in seq_along(collected)) {
            rownames(collected[[i]]) <- rownames(x)
            colnames(collected[[i]]) <- levels(f)
        }

    } else {
        lost <- is.na(ids)
        ids <- ids[!lost]
        if (any(lost)) {
            x <- x[,!lost,drop=FALSE]
        }

        # Drop unused levels, as the subsequent mapping step to preserve the type of 'ids'
        # in .create_coldata doesn't make sense (as there is no mapping to a concrete observation).
        by.group <- split(seq_along(ids), ids, drop=TRUE)

        out <- rowBlockApply(x, FUN=.summarize_assay_internal, by.group=by.group, 
            statistics=statistics, threshold=threshold, BPPARAM=BPPARAM)

        collected <- do.call(mapply, c(list(FUN=rbind, SIMPLIFY=FALSE, USE.NAMES=FALSE), out))
        names(collected) <- names(out[[1]])

        freq <- lengths(by.group)
        if ("mean" %in% statistics) {
            collected$mean <- t(t(collected$sum)/freq)
        }
        if ("prop.detected" %in% statistics) {
            collected$prop.detected <- t(t(collected$num.detected)/freq)
        }
    }

    list(summary=collected[statistics], freq=freq)
}

#' @importFrom MatrixGenerics rowSums rowMedians
#' @importClassesFrom SparseArray COO_SparseMatrix SVT_SparseMatrix
.summarize_assay_internal <- function(x, by.group, statistics, threshold) {
    if (is(x, "COO_SparseMatrix")) {
        x <- as(x, "SVT_SparseMatrix")
    }

    collated <- list()
    
    if ("sum" %in% statistics || "mean" %in% statistics) {
        out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE]))
        collated$sum <- .cbind_empty(out, x)
    }

    if ("median" %in% statistics) {
        out <- lapply(by.group, function(i) rowMedians(x[,i,drop=FALSE]))
        out <- .cbind_empty(out, x)
        rownames(out) <- rownames(x)
        collated$median <- out
    }

    if ("num.detected" %in% statistics || "prop.detected" %in% statistics) {
        out <- lapply(by.group, function(i) rowSums(x[,i,drop=FALSE] > threshold))
        collated$num.detected <- .cbind_empty(out, x)
    }

    collated
}

.cbind_empty <- function(out, x) {
    if (length(out)) {
        do.call(cbind, out)
    } else {
        as.matrix(x[,0,drop=FALSE])
    }
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
        coldata[[store.number]] <- unname(freq)
    }
    coldata
}

##########################
##########################

#' @export
#' @rdname summarizeAssayByGroup
setGeneric("summarizeAssayByGroup", function(x, ...) standardGeneric("summarizeAssayByGroup"))

#' @export
#' @rdname summarizeAssayByGroup
#' @importFrom BiocParallel SerialParam
setMethod("summarizeAssayByGroup", "ANY", .summarize_assay_by_group)

#' @export
#' @rdname summarizeAssayByGroup
#' @importFrom SummarizedExperiment assay
setMethod("summarizeAssayByGroup", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .summarize_assay_by_group(assay(x, assay.type), ...)
})
