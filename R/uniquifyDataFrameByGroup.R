#' Groupwise unique rows of a DataFrame
#'
#' Obtain unique values for groups of rows in a \linkS4class{DataFrame}.
#' This is used by \code{\link{aggregateAcrossCells}} to obtain \code{\link{colData}} for the aggregated SummarizedExperiment.
#'
#' @param x A \linkS4class{DFrame}.
#' @param grouping A factor (or a vector coercible into a factor) of length equal to \code{nrow(x)},
#' containing the groupings for each row of \code{x}.
#'
#' @return A \linkS4class{DFrame} where each row corresponds to a level of \code{grouping}.
#' Each column corresponds to a column of \code{x} and contains the unique value for each group, or \code{NA} if no such value exists.
#'
#' @author Aaron Lun
#' @examples
#' x <- DataFrame(
#'    foo = rep(LETTERS[1:3], 10),
#'    bar = rep(letters[1:3], 10),
#'    val = sample(3, 30, replace=TRUE)
#' )
#' uniquifyDataFrameByGroup(x, x$foo)
#' uniquifyDataFrameByGroup(x, x$val)
#'
#' @export
#' @importFrom S4Vectors DataFrame
uniquifyDataFrameByGroup <- function(x, grouping) {
    grouping <- factor(grouping)
    new.cd <- .merge_DF_rows(x, grouping, levels(grouping))
    do.call(DataFrame, c(new.cd, list(check.names=FALSE, row.names=levels(grouping))))
}
