#' Identify outlier values 
#' 
#' Convenience function to determine which values in a numeric vector are outliers based on the median absolute deviation (MAD).
#'
#' @param metric Numeric vector of values.
#' @param nmads A numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier.
#' @param type String indicating whether outliers should be looked for at both tails (\code{"both"}), only at the lower tail (\code{"lower"}) or the upper tail (\code{"higher"}).
#' @param log Logical scalar, should the values of the metric be transformed to the log2 scale before computing MADs?
#' @param subset Logical or integer vector, which subset of values should be used to calculate the median/MAD? 
#' If \code{NULL}, all values are used.
#' @param batch Factor of length equal to \code{metric}, specifying the batch to which each observation belongs. 
#' A median/MAD is calculated for each batch, and outliers are then identified within each batch.
#' @param share.medians Logical scalar indicating whether the median calculation should be shared across batches.
#' Only used if \code{batch} is specified.
#' @param share.mads Logical scalar indicating whether the MAD calculation should be shared across batches.
#' Only used if \code{batch} is specified.
#' @param share.missing Logical scalar indicating whether a common MAD/median should be used 
#' for any batch that has no values left after subsetting.
#' Only relevant when both \code{batch} and \code{subset} are specified.
#' @param min.diff A numeric scalar indicating the minimum difference from the median to consider as an outlier. 
#' Ignored if \code{NA}.
#' @param share_medians,share_mads,share_missing,min_diff
#' Soft-deprecated equivalents of the arguments above.
#' 
#' @return An outlier.filter object of the same length as the \code{metric} argument.
#' This is effectively a logical vector specifying the observations that are considered as outliers.
#' The chosen thresholds are stored in the \code{"thresholds"} attribute.
#'
#' @details
#' Lower and upper thresholds are stored in the \code{"thresholds"} attribute of the returned vector.
#' By default, this is a numeric vector of length 2 for the threshold on each side.
#' If \code{type="lower"}, the higher limit is \code{Inf}, while if \code{type="higher"}, the lower limit is \code{-Inf}.
#' 
#' If \code{min.diff} is not \code{NA}, the minimum distance from the median required to define an outlier is set as the larger of \code{nmads} MADs and \code{min.diff}.
#' This aims to avoid calling many outliers when the MAD is very small, e.g., due to discreteness of the metric.
#' If \code{log=TRUE}, this difference is defined on the log2 scale.
#' 
#' If \code{subset} is specified, the median and MAD are computed from a subset of cells and the values are used to define the outlier threshold that is applied to all cells.
#' In a quality control context, this can be handy for excluding groups of cells that are known to be low quality (e.g., failed plates) so that they do not distort the outlier definitions for the rest of the dataset.
#' 
#' Missing values trigger a warning and are automatically ignored during estimation of the median and MAD.
#' The corresponding entries of the output vector are also set to \code{NA} values.
#'
#' The outlier.filter class is derived from an ordinary logical vector.
#' The only difference is that any subsetting will not discard the \code{"thresholds"}, which avoids unnecessary loss of information.
#' Users can simply call \code{\link{as.logical}} to convert this into a logical vector.
#'
#' @section Handling batches:
#' If \code{batch} is specified, outliers are defined within each batch separately using batch-specific median and MAD values.
#' This gives the same results as if the input metrics were subsetted by batch and \code{isOutlier} was run on each subset,
#' and is often useful when batches are known \emph{a priori} to have technical differences (e.g., in sequencing depth).
#' 
#' If \code{share.medians=TRUE}, a shared median is computed across all cells.
#' If \code{share.mads=TRUE}, a shared MAD is computed using all cells 
#' (based on either a batch-specific or shared median, depending on \code{share.medians}).
#' These settings are useful to enforce a common location or spread across batches, e.g., we might set \code{share.mads=TRUE} for log-library sizes if coverage varies across batches but the variance across cells is expected to be consistent across batches.
#'
#' If a batch does not have sufficient cells to compute the median or MAD (e.g., after applying \code{subset}),
#' the default setting of \code{share.missing=TRUE} will set these values to the shared median and MAD.
#' This allows us to define thresholds for low-quality batches based on information in the rest of the dataset.
#' (Note that the use of shared values only affects this batch and not others unless \code{share.medians} and \code{share.mads} are also set.)
#' Otherwise, if \code{share.missing=FALSE}, all cells in that batch will have \code{NA} in the output.
#' 
#' If \code{batch} is specified, the \code{"threshold"} attribute in the returned vector is a matrix with one named column per level of \code{batch} and two rows (one per threshold).
#' 
#' @author Aaron Lun
#'
#' @examples
#' example_sce <- mockSCE()
#' stats <- perCellQCMetrics(example_sce)
#'
#' str(isOutlier(stats$sum))
#' str(isOutlier(stats$sum, type="lower"))
#' str(isOutlier(stats$sum, type="higher"))
#' 
#' str(isOutlier(stats$sum, log=TRUE))
#'
#' b <- sample(LETTERS[1:3], ncol(example_sce), replace=TRUE)
#' str(isOutlier(stats$sum, log=TRUE, batch=b))
#' 
#' @seealso
#' \code{\link{quickPerCellQC}}, a convenience wrapper to perform outlier-based quality control.
#'
#' \code{\link{perCellQCMetrics}}, to compute potential QC metrics.
#'
#' @aliases
#' outlier.filter
#' outlier.filter-class
#' [.outlier.filter
#' @export
isOutlier <- function(metric, nmads = 3, type = c("both", "lower", "higher"), 
    log = FALSE, subset = NULL, batch = NULL, share.medians=FALSE, 
    share.mads=FALSE, share.missing=TRUE, min.diff = NA,
    share_medians=NULL, share_mads=NULL, share_missing=NULL, min_diff=NULL)
{
    min.diff <- .replace(min.diff, min_diff)
    share.medians <- .replace(share.medians, share_medians)
    share.mads <- .replace(share.mads, share_mads)
    share.missing <- .replace(share.missing, share_missing)

    if (log) {
        metric <- log2(metric)
    }

    N <- length(metric)
    if (nobatch <- is.null(batch)) {
        batch <- rep("1", N)
    } else {
        if (length(batch) != N) { 
            stop("length of 'batch' must equal length of 'metric'")
        }
        batch <- as.character(batch)
    }

    stats <- .get_med_and_mad(metric, 
        batch=batch, 
        subset=subset, 
        share.medians=share.medians, 
        share.mads=share.mads,
        share.missing=share.missing)

    output <- .apply_filters(metric, 
        batch=batch, 
        cur.med=stats$med, 
        cur.mad=stats$mad,
        type=match.arg(type), 
        nmads=nmads, 
        min.diff=min.diff
    )

    thresholds <- output$thresholds
    if (nobatch) {
        thresholds <- drop(thresholds)
    }
    if (log) { 
        thresholds <- 2^thresholds
    }

    outliers <- outlier.filter(output$outliers)
    attr(outliers, "thresholds") <- thresholds
    outliers 
}

#' @importFrom stats mad median
.get_med_and_mad <- function(metric, batch, subset, share.medians, share.mads, share.missing) {
    # Subsetting by user-specific factor or by non-NA.
    if (!is.null(subset)) { 
        M <- metric[subset]
        B <- batch[subset]
    } else {
        M <- metric
        B <- batch
    }
    if (any(na.drop <- is.na(M))) { 
        M <- M[!na.drop]
        B <- B[!na.drop]
        warning("missing values ignored during outlier detection")
    }

    # Defining the mad or median, possibly by sharing information across batches.
    by.batch <- split(M, B)
    all_batches <- sort(unique(batch))
    empty <- rep(NA_real_, length(all_batches))
    names(empty) <- all_batches

    cur.med <- empty
    if (!share.medians) {
        cur.med[names(by.batch)] <- unlist(lapply(by.batch, median)) # handles no-batch input better than vapply.
    }
    replace <- .determine_sharingness(share.medians, share.missing, all_batches, names(by.batch))
    if (length(replace)) {
        cur.med[replace] <- median(M)
    }

    cur.mad <- empty
    if (!share.mads) {
        cur.mad[names(by.batch)] <- unlist(mapply(mad, x=by.batch, center=cur.med[names(by.batch)], SIMPLIFY=FALSE))
    }
    replace <- .determine_sharingness(share.mads, share.missing, all_batches, names(by.batch))
    if (length(replace)) {
        cur.mad[replace] <- median(abs(M - cur.med[B])) * formals(mad)$constant
    }

    list(med=cur.med, mad=cur.mad)
}

.apply_filters <- function(metric, batch, cur.med, cur.mad, type, nmads, min.diff) {
    diff.val <- pmax(min.diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val 
    lower.limit <- cur.med - diff.val 

    if (type == "lower") {
        upper.limit[] <- Inf
    } else if (type == "higher") {
        lower.limit[] <- -Inf
    }

    collected <- (metric < lower.limit[batch] | upper.limit[batch] < metric)
    names(collected) <- names(metric)

    all.threshold <- rbind(lower=lower.limit, higher=upper.limit)
    list(outliers=collected, thresholds=all.threshold)
}

.determine_sharingness <- function(share.value, share.missing, all_batches, valid_batches) {
    if (share.value && share.missing) {
        all_batches
    } else if (share.missing) {
        setdiff(all_batches, valid_batches)
    } else if (share.value) {
        valid_batches
    } else {
        character(0)
    }
}

#' @rawNamespace exportClasses(outlier.filter)
setOldClass(c("outlier.filter", "logical"))

#' @export
outlier.filter <- function(x) {
    class(x) <- c("outlier.filter", "logical")
    x
}

#' @export
#' @method [ outlier.filter
`[.outlier.filter` <- function(x, i, j, ..., drop=TRUE) {
    out <- NextMethod()
    at <- attributes(x)
    mostattributes(out) <- at[setdiff(names(at), "names")]
    out 
}
