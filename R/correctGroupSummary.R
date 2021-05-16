#' Correct group-level summaries
#'
#' Correct the summary statistic for each group for unwanted variation by fitting a linear model and extracting the coefficients.
#'
#' @param x A numeric matrix containing summary statistics for each gene (row) and combination of group and block (column),
#' computed by functions such as \code{\link{summarizeAssayByGroup}} - see Examples.
#' @param group A factor or vector specifying the group identity for each column of \code{x}, usually clusters or cell types.
#' @param block A factor or vector specifying the blocking level for each column of \code{x}, e.g., batch of origin.
#' @param transform String indicating how the differences between groups should be computed, for the batch adjustment.
#' @param offset Numeric scalar specifying the offset to use when \code{difference="log"} (default 1) or \code{difference="logit"} (default 0.01).
#' @param weights A numeric vector containing the weight of each combination, e.g., due to differences in the number of cells used to compute each summary.
#' If \code{NULL}, all combinations have equal weight.
#' 
#' @return A numeric matrix with number of rows equal to \code{nrow(x)} and number of columns equal to the number of unique levels in \code{group}.
#' Each column corresponds to a group and contains the averaged statistic across batches.
#'
#' @details
#' Here, we consider group-level summary statistics such as the average expression of all cells or the proportion with detectable expression.
#' These are easy to intepret and helpful for any visualizations that operate on individual groups, e.g., heatmaps.
#' 
#' However, in the presence of unwanted factors of variation (e.g., batch effects), some adjustment is required to ensure these group-level statistics are comparable.
#' We cannot directly average group-level statistics across batches as some groups may not exist in particular batches, e.g., due to the presence of unique cell types in different samples.
#' A direct average would be biased by variable contributions of the batch effect for each group.
#' 
#' To overcome this, we use groups that are present across multiple levels of the unwanted factor in multiple batches to correct for the batch effect.
#' (That is, any level of \code{groups} that occurs for multiple levels of \code{block}.)
#' For each gene, we fit a linear model to the (transformed) values containing both the group and block factors.
#' We then report the coefficient for each group as the batch-adjusted average for that group;
#' this is possible as the fitted model has no intercept.
#'
#' The default of \code{transform="raw"} will not transform the values, and is generally suitable for log-expression values.
#' Setting \code{transform="log"} will perform a log-transformation after adding \code{offset} (default of 1), and is suitable for normalized counts.
#' Setting \code{transform="logit"} will perform a logit transformation after adding \code{offset} (default of 0.01) - 
#' to the numerator and twice to the denominator, to shrink to 0.5 -
#' and is suitable for proportional data such as the proportion of detected cells.
#'
#' After the model is fitted to the transformed values, the reverse transformation is applied to the coefficients to obtain a corrected summary statistic on the original scale.
#' For \code{transform="log"}, any negative values are coerced to zero,
#' while for \code{transform="logit"}, any values outside of [0, 1] are coerced to the closest boundary.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{summarizeAssayByGroup}}, to generate the group-level summaries for this function.
#'
#' \code{regressBatches} from the \pkg{batchelor} package, to remove the batch effect from per-cell expression values.
#'
#' @examples
#' y <- matrix(rnorm(10000), ncol=1000)
#' group <- sample(10, ncol(y), replace=TRUE)
#' block <- sample(5, ncol(y), replace=TRUE)
#'
#' summaries <- summarizeAssayByGroup(y, DataFrame(group=group, block=block), 
#'     statistics=c("mean", "prop.detected"))
#'
#' # Computing batch-aware averages:
#' averaged <- correctGroupSummary(assay(summaries, "mean"), 
#'     group=summaries$group, block=summaries$block)
#' 
#' num <- correctGroupSummary(assay(summaries, "prop.detected"),
#'     group=summaries$group, block=summaries$block, transform="logit") 
#' 
#' @export
#' @importFrom stats lm.fit model.matrix lm.wfit
#' @importFrom Matrix t
correctGroupSummary <- function(x, group, block, transform=c("raw", "log", "logit"), offset=NULL, weights=NULL) { 
    transform <- match.arg(transform)
    if (transform=="log") {
        if (is.null(offset)) {
            offset <- 1
        }
        x <- log(x + offset)
    } else if (transform=="logit") {
        if (is.null(offset)) {
            offset <- 0.01
        }
        x <- (x + offset) / (1 + 2 * offset)
        x <- log(x/(1-x))
    }

    group <- factor(group)
    if (is.null(block) || length(unique(block))==1L) {
        design <- model.matrix(~0 + group)
    } else {
        design <- model.matrix(~0 + group + factor(block))
    }

    # Replace with fit.
    if (is.null(weights)) {
        fit <- lm.fit(y=t(x), x=design)
    } else {
        fit <- lm.wfit(y=t(x), x=design, w=weights)
    }

    if (nrow(x) > 1) {
        averages <- t(fit$coefficients[seq_len(nlevels(group)),,drop=FALSE])
    } else {
        averages <- rbind(head(fit$coefficients, nlevels(group)))
    }
    colnames(averages) <- levels(group)
    rownames(averages) <- rownames(x)

    if (transform=="log") {
        averages <- exp(averages) - offset
        averages[averages < 0] <- 0
    } else if (transform=="logit") {
        averages <- exp(averages)
        averages <- averages/(averages + 1)
        averages <- averages * (1 + 2 * offset) - offset
        averages[averages < 0] <- 0
        averages[averages > 1] <- 1
    }

    averages
}
