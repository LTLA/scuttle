#' Per-cell quality control metrics
#'
#' Compute per-cell quality control metrics for a count matrix or a \linkS4class{SingleCellExperiment}.
#'
#' @param x A numeric matrix of counts with cells in columns and features in rows.
#' 
#' Alternatively, a \linkS4class{SummarizedExperiment} or \linkS4class{SingleCellExperiment} object containing such a matrix.
#' @param subsets A named list containing one or more vectors 
#' (a character vector of feature names, a logical vector, or a numeric vector of indices),
#' used to identify interesting feature subsets such as ERCC spike-in transcripts or mitochondrial genes. 
#' @param percent.top An integer vector specifying the size(s) of the top set of high-abundance genes.
#' Used to compute the percentage of library size occupied by the most highly expressed genes in each cell.
#' @param threshold A numeric scalar specifying the threshold above which a gene is considered to be detected.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' @param ... For the generic, further arguments to pass to specific methods.
#' 
#' For the SummarizedExperiment and SingleCellExperiment methods, further arguments to pass to the ANY method.
#' @param assay.type A string or integer scalar indicating which \code{assays} in the \code{x} contains the count matrix.
#' @param use.altexps Logical scalar indicating whether QC statistics should be computed for alternative Experiments in \code{x}.
#' If \code{TRUE}, statistics are computed for all alternative experiments. 
#'
#' Alternatively, an integer or character vector specifying the alternative Experiments to use to compute QC statistics.
#' 
#' Alternatively \code{NULL}, in which case alternative experiments are not used.
#' @param flatten Logical scalar indicating whether the nested \linkS4class{DataFrame}s in the output should be flattened.
#' @param percent_top,detection_limit,exprs_values,use_altexps Soft deprecated equivalents to the arguments described above.
#'
#' @return
#' A \linkS4class{DataFrame} of QC statistics where each row corresponds to a column in \code{x}.
#' This contains the following fields:
#' \itemize{
#' \item \code{sum}: numeric, the sum of counts for each cell.
#' \item \code{detected}: numeric, the number of observations above \code{threshold}.
#' }
#'
#' If \code{flatten=FALSE}, the DataFrame will contain the additional columns:
#' \itemize{
#' \item \code{percent.top}: numeric matrix, the percentage of counts assigned to the top most highly expressed genes.
#' Each column of the matrix corresponds to an entry of \code{percent.top}, sorted in increasing order.
#' \item \code{subsets}: A nested DataFrame containing statistics for each subset, see Details.
#' \item \code{altexps}: A nested DataFrame containing statistics for each alternative experiment, see Details.
#' This is only returned for the SingleCellExperiment method.
#' \item \code{total}: numeric, the total sum of counts for each cell across main and alternative Experiments.
#' This is only returned for the SingleCellExperiment method.
#' }
#'
#' If \code{flatten=TRUE}, nested matrices and DataFrames are flattened to remove the hierarchical structure from the output DataFrame.
#' 
#' @author Aaron Lun
#' 
#' @details
#' This function calculates useful QC metrics for identification and removal of potentially problematic cells.
#' Obvious per-cell metrics are the sum of counts (i.e., the library size) and the number of detected features.
#' The percentage of counts in the top features also provides a measure of library complexity.
#' 
#' If \code{subsets} is specified, these statistics are also computed for each subset of features.
#' This is useful for investigating gene sets of interest, e.g., mitochondrial genes, Y chromosome genes.
#' These statistics are stored as nested \linkS4class{DataFrame}s in the \code{subsets} field of the output.
#' For example, if the input \code{subsets} contained \code{"Mito"} and \code{"Sex"}, the output would look like:
#' \preformatted{  output 
#'   |-- sum
#'   |-- detected
#'   |-- percent.top
#'   +-- subsets
#'       |-- Mito
#'       |   |-- sum
#'       |   |-- detected
#'       |   +-- percent
#'       +-- Sex 
#'           |-- sum
#'           |-- detected
#'           +-- percent
#' }
#' Here, the \code{percent} field contains the percentage of each cell's count sum assigned to each subset. 
#'
#' If \code{use.altexps=TRUE}, the same statistics are computed for each alternative experiment in \code{x}.
#' This can also be an integer or character vector specifying the alternative Experiments to use.
#' These statistics are also stored as nested \linkS4class{DataFrame}s, this time in the \code{altexps} field of the output.
#' For example, if \code{x} contained the alternative Experiments \code{"Spike"} and \code{"Ab"}, the output would look like:
#' \preformatted{  output 
#'   |-- sum
#'   |-- detected
#'   |-- percent.top
#'   +-- altexps 
#'   |   |-- Spike
#'   |   |   |-- sum
#'   |   |   |-- detected
#'   |   |   +-- percent.total
#'   |   +-- Ab
#'   |       |-- sum
#'   |       |-- detected
#'   |       +-- percent.total
#'   +-- total 
#' }
#' The \code{total} field contains the total sum of counts for each cell across the main and alternative Experiments.
#' The \code{percent} field contains the percentage of the total count in each alternative Experiment for each cell.
#' 
#' If \code{flatten=TRUE}, the nested DataFrames are flattened by concatenating the column names with underscores.
#' This means that, say, the \code{subsets$Mito$sum} nested field becomes the top-level \code{subsets_Mito_sum} field.
#' A flattened structure is more convenient for end-users performing interactive analyses,
#' but less convenient for programmatic access as artificial construction of strings is required.
#' 
#' @examples
#' example_sce <- mockSCE()
#' stats <- perCellQCMetrics(example_sce)
#' stats
#'
#' # With subsets.
#' stats2 <- perCellQCMetrics(example_sce, subsets=list(Mito=1:10), 
#'     flatten=FALSE)
#' stats2$subsets
#'
#' # With alternative Experiments.
#' pretend.spike <- ifelse(seq_len(nrow(example_sce)) < 10, "Spike", "Gene")
#' alt_sce <- splitAltExps(example_sce, pretend.spike)
#' stats3 <- perCellQCMetrics(alt_sce, flatten=FALSE)
#' stats3$altexps
#'
#'
#' @seealso 
#' \code{\link{addPerCellQC}}, to add the QC metrics to the column metadata.
#' @export
#' @name perCellQCMetrics
NULL

#' @importFrom beachmat colBlockApply
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
#' @importFrom BiocParallel bpmapply SerialParam
#' @importClassesFrom S4Vectors DataFrame
.per_cell_qc_metrics <- function(x, subsets = NULL, percent.top = integer(0),
    threshold = 0, BPPARAM=SerialParam(), flatten=TRUE, 
    percent_top=NULL, detection_limit=NULL)
{
    threshold <- .replace(threshold, detection_limit)
    percent.top <- .replace(percent.top, percent_top)

    if (length(subsets) && is.null(names(subsets))){ 
        stop("'subsets' must be named")
    }
    subsets <- lapply(subsets, FUN = .subset2index, target = x, byrow = TRUE)
    percent.top <- sort(as.integer(percent.top))

    # Computing all QC metrics, with cells split across workers. 
    bp.out <- colBlockApply(x, FUN=.per_cell_qc, featcon=subsets, top=percent.top, limit=threshold, BPPARAM=BPPARAM)

    # Aggregating across cores.
    full.info <- DataFrame(
        sum=unlist(lapply(bp.out, FUN=function(x) x[[1]][[1]])),
        detected=unlist(lapply(bp.out, FUN=function(x) x[[1]][[2]])),
        row.names=colnames(x)
    )

    pct <- do.call(cbind, lapply(bp.out, FUN=function(x) x[[1]][[3]]))
    rownames(pct) <- percent.top
    full.info$percent.top <- t(pct)/full.info$sum * 100

    # Collecting subset information.
    if (!is.null(subsets)) {
        sub.info <- make_zero_col_DFrame(ncol(x))
        for (i in seq_along(subsets)) {
            sub.out <- DataFrame(
                sum=unlist(lapply(bp.out, FUN=function(x) x[[2]][[i]][[1]])),
                detected=unlist(lapply(bp.out, FUN=function(x) x[[2]][[i]][[2]]))
            )
            sub.out$percent <- sub.out$sum/full.info$sum * 100
            sub.info[[names(subsets)[i]]] <- sub.out
        }
        full.info$subsets <- sub.info
    }

    if (flatten) {
        full.info <- .flatten_nested_dims(full.info)
    }
    full.info
}

#' @importFrom Matrix colSums 
#' @importClassesFrom Matrix sparseMatrix
#' @importClassesFrom DelayedArray SparseArraySeed
.per_cell_qc <- function(x, featcon, top, limit) {
    if (is(x, "SparseArraySeed")) {
        x <- as(x, "sparseMatrix")
    }

    detected <- x > limit

    full <- list(
        sum=unname(colSums(x)),
        detected=unname(colSums(detected)),
        prop=cumulative_prop(x, top)
    )

    featcons <- lapply(featcon, function(i) {
        # TODO: switch to MatrixGenerics when that finally becomes available.
        list(
            sum=unname(colSums(x[i,,drop=FALSE])),
            detected=unname(colSums(detected[i,,drop=FALSE]))
        )
    })

    list(full, featcons)
}

##################################################

#' @export
#' @rdname perCellQCMetrics
setGeneric("perCellQCMetrics", function(x, ...) standardGeneric("perCellQCMetrics"))

#' @export
#' @rdname perCellQCMetrics
setMethod("perCellQCMetrics", "ANY", .per_cell_qc_metrics)

#' @export
#' @rdname perCellQCMetrics
#' @importFrom SummarizedExperiment assay
setMethod("perCellQCMetrics", "SummarizedExperiment", function(x, ..., assay.type="counts", exprs_values=NULL) {
    assay.type <- .replace(assay.type, exprs_values)
    .per_cell_qc_metrics(assay(x, assay.type), ...)
})

#' @export
#' @rdname perCellQCMetrics
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp altExpNames
#' @importFrom S4Vectors make_zero_col_DFrame I
#' @importClassesFrom S4Vectors DataFrame
setMethod("perCellQCMetrics", "SingleCellExperiment", 
    function(x, subsets=NULL, percent.top=integer(0), ..., flatten=TRUE, assay.type="counts", use.altexps=TRUE, 
        percent_top=NULL, exprs_values=NULL, use_altexps=NULL) 
{
    assay.type <- .replace(assay.type, exprs_values)
    use.altexps <- .replace(use.altexps, use_altexps)
    percent.top <- .replace(percent.top, percent_top)

    # subsets and percent.top need to be explicitly listed,
    # because the altexps call sets them to NULL and integer(0).
    main <- .per_cell_qc_metrics(assay(x, assay.type), subsets=subsets, percent.top=percent.top, flatten=FALSE, ...)
    use.altexps <- .use_names_to_integer_indices(use.altexps, x=x, nameFUN=altExpNames, msg="use.altexps")

    alt <- list()
    total <- main$sum
    for (i in seq_along(use.altexps)) {
        y <- assay(altExp(x, use.altexps[i]), assay.type)
        current <- .per_cell_qc_metrics(y, subsets=NULL, percent.top=integer(0), ...)
        current$percent.top <- current$subsets <- NULL
        total <- total + current$sum
        alt[[i]] <- current
    }
    for (i in seq_along(alt)) {
        alt[[i]]$percent <- alt[[i]]$sum/total * 100
    }

    if (length(alt)) {
        main$altexps <- do.call(DataFrame, lapply(alt, I))
        names(main$altexps) <- altExpNames(x)[use.altexps]
    } else {
        main$altexps <- make_zero_col_DFrame(ncol(x))
    }

    main$total <- total
    if (flatten) {
        main <- .flatten_nested_dims(main)
    }
    main
})

##################################################

#' @importFrom S4Vectors DataFrame
.flatten_nested_dims <- function(x, name="") {
    if (!is.null(dim(x))) {
        if (name!="") {
            name <- paste0(name, "_")
        }
        names <- sprintf("%s%s", name, colnames(x))
        rn <- rownames(x)

        df <- vector("list", ncol(x))
        for (i in seq_along(df)) {
            df[[i]] <- .flatten_nested_dims(x[,i], names[i])
        }
        if (length(df) > 0) {
            df <- do.call(cbind, df)
        } else {
            df <- DataFrame(x[,0])
        }

        rownames(df) <- rn
    } else {
        df <- DataFrame(x)
        colnames(df) <- name
    }
    df
}
