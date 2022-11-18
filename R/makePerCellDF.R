#' Create a per-cell data.frame
#'
#' Create a per-cell data.frame (i.e., where each row represents a cell) from a \linkS4class{SingleCellExperiment},
#' most typically for creating custom \pkg{ggplot2} plots.
#'
#' @param x A \linkS4class{SingleCellExperiment} object.
#' This is expected to have non-\code{NULL} row names.
#' @param features Character vector specifying the features for which to extract expression profiles across cells.
#' May also include features in alternative Experiments if permitted by \code{use_altexps}.
#' @param assay.type String or integer scalar indicating the assay to use to obtain expression values.
#' Must refer to a matrix-like object with integer or numeric values.
#' @param use.coldata Logical scalar indicating whether column metadata of \code{x} should be included.
#' Alternatively, a character or integer vector specifying the column metadata fields to use.
#' @param use.altexps Logical scalar indicating whether (meta)data should be extracted for alternative experiments in \code{x}.
#' Alternatively, a character or integer vector specifying the alternative experiments to use.
#' @param use.dimred Logical scalar indicating whether data should be extracted for dimensionality reduction results in \code{x}.
#' Alternatively, a character or integer vector specifying the dimensionality reduction results to use.
#' @param prefix.altexps Logical scalar indicating whether \code{\link{altExp}}-derived fields should be prefixed with the name of the alternative Experiment.
#' @param check.names Logical scalar indicating whether column names of the output should be made syntactically valid and unique.
#' @param exprs_values,use_dimred,use_altexps,prefix_altexps,check_names
#' Soft-deprecated equivalents of the arguments described above.
#'
#' @return A data.frame containing one field per aspect of data in \code{x} - see Details.
#' Each row corresponds to a cell (i.e., column) of \code{x}.
#'
#' @details
#' This function enables us to conveniently create a per-feature data.frame from a \linkS4class{SingleCellExperiment}.
#' Each row of the returned data.frame corresponds to a column in \code{x},
#' while each column of the data.frame corresponds to one aspect of the (meta)data in \code{x}.
#'
#' Columns are provided in the following order:
#' \enumerate{
#' \item Columns named according to the entries of \code{features} represent the expression values across cells for the specified feature in the \code{assay.type} assay.
#' \item Columns named according to the columns of \code{colData(x)} represent column metadata variables.
#' This consists of all variables if \code{use.coldata=TRUE}, no variables if \code{use.coldata=FALSE},
#' and only the specified variables if \code{use.coldata} is set to an integer or character vector.
#' \item Columns named in the format of \code{<DIM>.<NUM>} represent the \code{<NUM>}th dimension of the dimensionality reduction result \code{<DIM>}.
#' This is generated for all dimensionality reduction results if \code{use.dimred=TRUE}, none if \code{use.dimred=FALSE},
#' and only the specified results if \code{use.dimred}is set to an integer or character vector.
#' \item Columns named according to the row names of successive alternative Experiments,
#' representing the assay data in these objects.
#' These columns are only included if they are specified in \code{features} and if \code{use.altexps} is set.
#' Column names are prefixed with the name of the alternative Experiment if \code{prefix.altexps=TRUE}.
#' }
#'
#' By default, nothing is done to resolve syntactically invalid or duplicated column names.
#' \code{check_names=TRUE}, this is resolved by passing the column names through \code{\link{make.names}}.
#' Of course, as a result, some columns may not have the same names as the original fields in \code{x}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{makePerFeatureDF}}, for the feature-level equivalent.
#'
#' @examples
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' reducedDim(sce, "PCA") <- matrix(rnorm(ncol(sce)*10), ncol=10) # made-up PCA.
#'
#' df <- makePerCellDF(sce, features="Gene_0001")
#' head(df)
#'
#' @export
#' @importFrom SingleCellExperiment colData reducedDims reducedDimNames altExps altExpNames
makePerCellDF <- function(x, features=NULL, assay.type="logcounts",
    use.coldata=TRUE, use.dimred=TRUE, use.altexps=TRUE, prefix.altexps=FALSE, check.names=FALSE,
    swap_rownames = NULL, exprs_values=NULL, use_dimred=NULL, use_altexps=NULL, prefix_altexps=NULL,
    check_names=NULL)
{
    use.dimred <- .replace(use.dimred, use_dimred)
    use.altexps <- .replace(use.altexps, use_altexps)
    assay.type <- .replace(assay.type, exprs_values)
    prefix.altexps <- .replace(prefix.altexps, prefix_altexps)
    check.names <- .replace(check.names, check_names)

    # Initialize output list
    output <- list()

    # Collecting the column metadata.
    use.coldata <- .use_names_to_integer_indices(use.coldata, x=x, nameFUN=function(x) colnames(colData(x)), msg="use.coldata")
    if (length(use.coldata)) {
        cd <- colData(x)[,use.coldata,drop=FALSE]
        output <- c(output, list(as.data.frame(cd)))
    }

    # Collecting feature data
    output <- c(output, .harvest_se_by_column(x, features=features, assay.type=assay.type, swap_rownames = swap_rownames))

    # Collecting the reduced dimensions.
    use.dimred <- .use_names_to_integer_indices(use.dimred, x=x, nameFUN=reducedDimNames, msg="use.dimred")
    if (length(use.dimred)) {
        all_reds <- reducedDims(x)[use.dimred]
        red_vals <- vector("list", length(all_reds))

        for (r in seq_along(red_vals)) {
            curred <- data.frame(all_reds[[r]])
            names(curred) <- sprintf("%s.%s", names(all_reds)[r], seq_len(ncol(curred)))
            red_vals[[r]] <- curred
        }

        red_vals <- do.call(cbind, red_vals)
        output <- c(output, list(red_vals))
    }

    # Collecting the alternative Experiments.
    use.altexps <- .use_names_to_integer_indices(use.altexps, x=x, nameFUN=altExpNames, msg="use.altexps")
    if (length(use.altexps)) {
        all_alts <- altExps(x)[use.altexps]
        alt_vals <- vector("list", length(all_alts))

        for (a in seq_along(alt_vals)) {
            curalt <- .harvest_se_by_column(all_alts[[a]], features=features, assay.type=assay.type, swap_rownames = swap_rownames)
            if (prefix.altexps) {
                colnames(curalt) <- sprintf("%s.%s", names(all_alts)[a], colnames(curalt))
            }
            alt_vals[[a]] <- curalt
        }

        alt_vals <- do.call(cbind, alt_vals)
        output <- c(output, list(alt_vals))
    }

    # Checking the names.
    output <- do.call(cbind, output)
    if (check.names) {
        colnames(output) <- make.names(colnames(output), unique=TRUE)
    }
    output
}

#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix t
.harvest_se_by_column <- function(x, features, assay.type, swap_rownames = NULL) {
    x <- .swap_rownames(x, swap_rownames = swap_rownames)
    keep <- rownames(x) %in% features
    if (any(keep)) {
        curmat <- assay(x, assay.type, withDimnames=FALSE)[keep,,drop=FALSE]
        curmat <- as.matrix(t(curmat))
        assay_vals <- data.frame(curmat, row.names=colnames(x))
        colnames(assay_vals) <- rownames(x)[keep]
        assay_vals
    } else {
        # Avoid throwing an error for altexps if the feature doesn't even match.
        data.frame(matrix(0, ncol(x), 0L), row.names=colnames(x))
    }
}


## Using non-standard gene IDs
# from scater
.swap_rownames <- function(object, swap_rownames = NULL) {
    if (is.null(swap_rownames)) {
        return(object)
    }
    rownames(object) <- .get_rowData_column(object, swap_rownames)
    object
}

# from scater
.get_rowData_column <- function(object, column) {
    if (!column %in% colnames(rowData(object))) {
        stop("Cannot find column ", column, " in rowData")
    }
    rowData(object)[[column]]
}

