#' Normalization by deconvolution
#'
#' Scaling normalization of single-cell RNA-seq data by deconvolving size factors from cell pools.
#' 
#' @param x For \code{pooledSizeFactors}, a  numeric matrix-like object of counts, where rows are genes and columns are cells.
#' Alternatively, a \linkS4class{SummarizedExperiment} object containing such a matrix.
#'
#' For \code{computePooledFactors}, a \linkS4class{SingleCellExperiment} object containing a count matrix.
#' @param sizes A numeric vector of pool sizes, i.e., number of cells per pool.
#' @param clusters An optional factor specifying which cells belong to which cluster, for deconvolution within clusters.
#' @param ref.clust A level of \code{clusters} to be used as the reference cluster for inter-cluster normalization.
#' @param max.cluster.size An integer scalar specifying the maximum number of cells in each cluster.
#' @param positive A logical scalar indicating whether linear inverse models should be used to enforce positive estimates.
#' @param scaling A numeric scalar containing scaling factors to adjust the counts prior to computing size factors.
#' @param min.mean A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.
#' @param subset.row An integer, logical or character vector specifying the features to use.
#' @param BPPARAM A BiocParallelParam object specifying whether and how clusters should be processed in parallel.
#' @param ... For the \code{pooledSizeFactors} generic, additional arguments to pass to each method.
#' For the \linkS4class{SummarizedExperiment} method, additional methods to pass to the ANY method.
#' 
#' For the \code{computePooledFactors} function, additional arguments to pass to \code{pooledSizeFactors}.
#' @param assay.type A string specifying which assay values to use when \code{x} is a SummarizedExperiment or SingleCellExperiment.
#' 
#' @section Overview of the deconvolution method:
#' The \code{pooledSizeFactors} function implements the deconvolution strategy of Lun et al. (2016) for scaling normalization of sparse count data.
#' Briefly, a pool of cells is selected and the expression profiles for those cells are summed together.
#' The pooled expression profile is normalized against an average reference pseudo-cell, constructed by averaging the counts across all cells.
#' This defines a size factor for the pool as the median ratio between the count sums and the average across all genes.
#' 
#' The scaling bias for the pool is equal to the sum of the biases for the constituent cells.
#' The same applies for the size factors, as these are effectively estimates of the bias for each cell.
#' This means that the size factor for the pool can be written as a linear equation of the size factors for the cells.
#' Repeating this process for multiple pools will yield a linear system that can be solved to obtain the size factors for the individual cells.
#' 
#' In this manner, pool-based factors are deconvolved to yield the relevant cell-based factors.
#' The advantage is that the pool-based estimates are more accurate, as summation reduces the number of stochastic zeroes and the associated bias of the size factor estimate.
#' This accuracy feeds  back into the deconvolution process, thus improving the accuracy of the cell-based size factors.
#' 
#' @section Pooling with a sliding window:
#' Within each cluster (if not specified, all cells are put into a single cluster), cells are sorted by increasing library size and a sliding window is applied to this ordering.
#' Each location of the window defines a pool of cells with similar library sizes.
#' This avoids inflated estimation errors for very small cells when they are pooled with very large cells.
#' Sliding the window will construct an over-determined linear system that can be solved by least-squares methods to obtain cell-specific size factors.
#' 
#' Window sliding is repeated with different window sizes to construct the linear system, as specified by \code{sizes}.
#' By default, the number of cells in each window ranges from 21 to 101.
#' Using a range of window sizes improves the precision of the estimates, at the cost of increased computational work.
#' The defaults were chosen to provide a reasonable compromise between these two considerations.
#' The default set of \code{sizes} also avoids rare cases of linear dependencies and unstable estimates when all pool sizes are not co-prime with the number of cells.
#' 
#' The smallest window should be large enough so that the pool-based size factors are, on average, non-zero.
#' We recommend window sizes no lower than 20 for UMI data, though smaller windows may be possible for read count data.
#' The total number of cells should also be at least 100 for effective pooling.
#' (If \code{cluster} is specified, we would want at least 100 cells per cluster.)
#' 
#' If there are fewer cells than the smallest window size, the function will naturally degrade to performing library size normalization.
#' This yields results that are the same as \code{\link{librarySizeFactors}}.
#' 
#' @section Prescaling of the counts:
#' The simplest approach to pooling is to simply add the counts together for all cells in each pool.
#' However, this is suboptimal as any errors in the estimation of the pooled size factor will propagate to all component cell-specific size factors upon solving the linear system.
#' If the error is distributed evenly across all cell-specific size factors, the small size factors will have larger relative errors compared to the large size factors.
#' 
#' To avoid this, we perform \dQuote{prescaling} where we divide the counts by a cell-specific factor prior to pooling.
#' Ideally, the prescaling factor should be close to the true size factor for each cell.
#' Solving the linear system constructed with prescaled values should yield estimates that are more-or-less equal across all cells.
#' Thus, given similar absolute errors, the relative errors for all cells will also be similar.
#'
#' Obviously, the true size factor is unknown (otherwise why bother running this function?)
#' so we use the library size for each cell as a proxy instead.
#' This may perform poorly in pathological scenarios involving extreme differential expression and strong composition biases.
#' In cases where a more appropriate initial estimate is available, 
#' this can be used as the prescaling factor by setting the \code{scaling} argument.
#'
#' One potential approach is to run \code{computePooledFactors} twice to improve accuracy.
#' The first run is done as usual and will yield an initial estimate of the size factor for each cell.
#' In the second run, we supply our initial estimates in the \code{scaling} argument to serve as better prescaling factors.
#' Obviously, this involves twice as much computational work so we would only recommend attempting this in extreme circumstances.
#' 
#' @section Solving the linear system:
#' The linear system is solved using the sparse QR decomposition from the \pkg{Matrix} package.
#' However, this has known problems when the linear system becomes too large (see \url{https://stat.ethz.ch/pipermail/r-help/2011-August/285855.html}).
#' In such cases, we set \code{clusters} to break up the linear system into smaller, more manageable components that can be solved separately.
#' The default \code{max.cluster.size} will arbitrarily break up the cell population (within each cluster, if specified) so that we never pool more than 3000 cells.
#' Note that this involves appending a suffix like \code{"-1"} to the end of each cluster's name;
#' this may appear on occasion in warnings or error messages.
#' 
#' @section Normalization within and between clusters:
#' In general, it is more appropriate to pool more similar cells to avoid violating the assumption of a non-DE majority of genes.
#' This can be done by specifying the \code{clusters} argument where cells in each cluster have similar expression profiles.
#' Deconvolution is subsequently applied on the cells within each cluster, where there should be fewer DE genes between cells.
#' Any clustering can be used, and only a rough clustering is required; \code{computePooledFactors} is robust to a moderate level of DE within each cluster.
#' The \code{\link[scran]{quickCluster}} function from the \pkg{scran} package is particularly convenient for this purpose.
#' 
#' Size factors computed within each cluster must be rescaled for comparison between clusters.
#' To do so, we choose one cluster as a \dQuote{reference} to which all others are normalized.
#' Ideally, the reference cluster should have a stable expression profile and not be extremely different from all other clusters.
#' The assumption here is that there is a non-DE majority between the reference and each other cluster
#' (which is still a weaker assumption than that required without clustering).
#' The rescaling factor is then defined by computing the ratios in averaged expression between each cluster's pseudo-cell and that of the reference,
#' and taking the median of these ratios across all genes.
#' 
#' By default, the cluster with the most non-zero counts is used as the reference.
#' This reduces the risk of obtaining undefined rescaling factors for the other clusters, while improving the precision (and also accuracy) of the median-based factor estimate.
#' Alternatively, the reference can be manually specified using \code{ref.clust} if there is prior knowledge about which cluster is most suitable, e.g., from PCA or t-SNE plots.
#'
#' Each cluster should ideally be large enough to contain a sufficient number of cells for pooling.
#' Otherwise, \code{computePooledFactors} will fall back to library size normalization for small clusters.
#' 
#' If the estimated rescaling factor is not positive, a warning is emitted and the function falls back to the ratio of sums between pseudo-cells (in effect, library size normalization).
#' This can occasionally happen when a cluster's cells expresses a small subset of genes - 
#' this is not problematic for within-cluster normalization, as non-expressed genes are simply ignored,
#' but violates the assumption of a non-DE majority when performing inter-cluster comparisons.
#' 
#' @section Dealing with non-positive size factors:
#' It is possible for the deconvolution algorithm to yield negative or zero estimates for the size factors.
#' These values are obviously nonsensical and \code{computePooledFactors} will raise a warning if they are encountered.
#' Negative estimates are mostly commonly generated from low quality cells with few expressed features, such that most genes still have zero counts even after pooling.
#' They may also occur if insufficient filtering of low-abundance genes was performed.
#' 
#' To avoid these problematic size factors, the best solution is to increase the stringency of the filtering.
#' \itemize{
#' \item If only a few negative/zero size factors are present, they are likely to correspond to a few low-quality cells with few expressed features.
#' Such cells are difficult to normalize reliably under any approach, and can be removed by increasing the stringency of the quality control.
#' \item If many negative/zero size factors are present, it is probably due to insufficient filtering of low-abundance genes.
#' This results in many zero counts and pooled size factors of zero, and can be fixed by filtering out more genes with a higher \code{min.mean} - see \dQuote{Gene selection} below.
#' }
#' Another approach is to increase in the number of \code{sizes} to improve the precision of the estimates.
#' This reduces the chance of obtaining negative/zero size factors due to estimation error, for cells where the true size factors are very small.
#' 
#' As a last resort, \code{positive=TRUE} is set by default, which uses \code{\link{cleanSizeFactors}} to coerce any non-positive estimates to positive values.
#' This ensures that, at the very least, downstream analysis is possible even if the size factors for affected cells are not accurate.
#' Users can skip this step by setting \code{positive=FALSE} to perform their own diagnostics or coercions.
#' 
#' @section Gene selection:
#' If too many genes have consistently low counts across all cells, even the pool-based size factors will be zero.
#' This results in zero or negative size factor estimates for many cells.
#' We avoid this by filtering out low-abundance genes using the \code{min.mean} argument.
#' This represents a minimum threshold \code{min.mean} on the library size-adjusted average counts from \code{\link{calculateAverage}}.
#' 
#' By default, we set \code{min.mean} to 1 for read count data and 0.1 for UMI data.
#' The exact values of these defaults are more-or-less arbitrary and are retained for historical reasons.
#' The lower threshold for UMIs is motivated by (i) their lower count sizes, which would result in the removal of too many genes with a higher threshold; and (ii) the lower variability of UMI counts, which results in a lower frequency of zeroes compared to read count data at the same mean.
#' We use the median library size to detect whether the counts are those of reads (above 100,000) or UMIs (below 50,000) to automatically set \code{min.mean}.
#' Mean library sizes in between these two limits will trigger a warning and revert to using \code{min.mean=0.1}.
#' 
#' If \code{clusters} is specified, filtering by \code{min.mean} is performed on the per-cluster average during within-cluster normalization,
#' and then on the (library size-adjusted) average of the per-cluster averages during between-cluster normalization.
#'
#' Performance can generally be improved by removing genes that are known to be strongly DE between cells.
#' This weakens the assumption of a non-DE majority and avoids the error associated with DE genes.
#' For example, we might remove viral genes when our population contains both infected and non-infected cells.
#' Of course, \code{computePooledFactors} is robust to some level of DE genes - that is, after all, its raison d'etre -
#' so one should only explicitly remove DE genes if it is convenient to do so. 
#' 
#' @section Obtaining standard errors:
#' Previous versions of \code{computePooledFactors} would return the standard error for each size factor when \code{errors=TRUE}.
#' This argument is no longer available as we have realized that standard error estimation from the linear model is not reliable.
#' Errors are likely underestimated due to correlations between pool-based size factors when they are computed from a shared set of underlying counts.
#' Users wishing to obtain a measure of uncertainty are advised to perform simulations instead, using the original size factor estimates to scale the mean counts for each cell.
#' Standard errors can then be calculated as the standard deviation of the size factor estimates across simulation iterations.
#' 
#' @return
#' For \code{pooledSizeFactors}, a numeric vector of size factors for all cells in \code{x} is returned.
#' 
#' For \code{computePooledFactors}, an object of class \code{x} is returned containing the vector of size factors in \code{\link{sizeFactors}(x)}.
#' 
#' @author
#' Aaron Lun and Karsten Bach
#' 
#' @seealso
#' \code{\link{logNormCounts}}, which uses the computed size factors to compute normalized expression values.
#'
#' \code{\link{librarySizeFactors}} and \code{\link{medianSizeFactors}}, for simpler approaches to computing size factors.
#'
#' \code{\link[scran]{quickCluster}} from the \pkg{scran} package, to obtain a rough clustering for use in \code{clusters}.
#'
#' @examples
#' library(scuttle)
#' sce <- mockSCE(ncells=500)
#' 
#' # Computing the size factors.
#' sce <- computePooledFactors(sce)
#' head(sizeFactors(sce))
#' plot(librarySizeFactors(sce), sizeFactors(sce), log="xy")
#'
#' # Using pre-clustering.
#' library(scran)
#' preclusters <- quickCluster(sce)
#' table(preclusters)
#' 
#' sce2 <- computePooledFactors(sce, clusters=preclusters)
#' head(sizeFactors(sce2))
#'
#' @references
#' Lun ATL, Bach K and Marioni JC (2016).
#' Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
#' \emph{Genome Biol.} 17:75
#'
#' @name computePooledFactors
NULL

#' @importFrom BiocParallel bplapply SerialParam
.calculate_pooled_factors <- function(x, sizes=seq(21, 101, 5), clusters=NULL, ref.clust=NULL, max.cluster.size=3000, 
    positive=TRUE, scaling=NULL, min.mean=NULL, subset.row=NULL, BPPARAM=SerialParam())
# This contains the function that performs normalization on the summed counts.
# It also provides support for normalization within clusters, and then between
# clusters to make things comparable. 
{
    ncells <- ncol(x)
    if (is.null(clusters)) {
        clusters <- integer(ncells)
    }
	clusters <- .limit_cluster_size(clusters, max.cluster.size)

    if (ncells!=length(clusters)) { 
        stop("'ncol(x)' is not equal to 'length(clusters)'")
    }
    indices <- split(seq_along(clusters), clusters)

    if (length(indices)==0L || any(lengths(indices)==0L)) {
        stop("zero cells in one of the clusters")
    }

    # Addigional sanity checks on various parameters.
    if (!is.null(scaling) && length(scaling)!=ncol(x)) {
        stop("'length(scaling)' should be equal to 'ncol(x)'")
    }

    min.mean <- .guessMinMean(x, min.mean=min.mean, BPPARAM=BPPARAM)

    sizes <- sort(as.integer(sizes))
    if (anyDuplicated(sizes)) { 
        stop("'sizes' are not unique") 
    }

    # Fragmenting the matrices (and also scaling).
    frag.x <- frag.scale <- vector("list", length(indices))
    for (i in seq_along(indices)) {
        idx <- indices[[i]]
        if (length(indices) > 1L || !identical(idx, seq_along(idx))) {
            current <- x[,idx,drop=FALSE]
        } else {
            current <- x
        }
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }
        frag.x[[i]] <- current
        frag.scale[i] <- list(scaling[idx]) # handle NULLs properly.
    }

    # Computing normalization factors within each cluster.
    all.norm <- bpmapply(FUN=.per_cluster_normalize, x=frag.x, scaling=frag.scale, 
        MoreArgs=list(sizes=sizes, min.mean=min.mean, positive=positive),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    names(all.norm) <- names(indices)

    clust.nf <- lapply(all.norm, "[[", i="final.nf")
    clust.profile <- lapply(all.norm, "[[", i="ave.cell")

    # Adjusting size factors between clusters.
    if (is.null(ref.clust)) {
        non.zeroes <- vapply(clust.profile, FUN=function(x) sum(x>0), FUN.VALUE=0L) 
        ref.clust <- which.max(non.zeroes)
    }
    rescaling.factors <- .rescale_clusters(clust.profile, ref.col=ref.clust, min.mean=min.mean) 

    clust.nf.scaled <- Map(`*`, clust.nf, rescaling.factors)
    clust.nf.scaled <- unlist(clust.nf.scaled)

    # Returning centered size factors, rather than normalization factors.
    final.sf <- rep(NA_real_, ncells)
    indices <- unlist(indices)
    final.sf[indices] <- clust.nf.scaled
    
    is.pos <- final.sf > 0 & !is.na(final.sf)
    final.sf/mean(final.sf[is.pos])
}

#' @export
#' @importFrom stats median
#' @importFrom MatrixGenerics colSums
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.guessMinMean <- function(x, min.mean, BPPARAM) { 
    # Choosing a mean filter based on the data type and then filtering:
    if (is.null(min.mean)) {
        old <- getAutoBPPARAM()
        setAutoBPPARAM(BPPARAM)
        on.exit(setAutoBPPARAM(old))

        mid.lib <- median(colSums(x))
        if (is.na(mid.lib)) { # no column check, for safety.
            min.mean <- 1
        } else if (mid.lib <= 50000) { # Probably UMI data.
            min.mean <- 0.1
        } else if (mid.lib >= 100000) { # Probably read data.
            min.mean <- 1
        } else {
            min.mean <- 0.1
            warning("assuming UMI data when setting 'min.mean'")
        }
    } else {
        min.mean <- pmax(min.mean, 1e-8) # must be positive.
    }
    min.mean
}

#############################################################
# Internal functions.
#############################################################

#' @importFrom Matrix qr qr.coef
#' @importFrom S4Arrays is_sparse
#' @importFrom MatrixGenerics colSums
.per_cluster_normalize <- function(x, sizes, min.mean=NULL, positive=FALSE, scaling=NULL) 
# Computes the normalization factors _within_ each cluster,
# along with the reference pseudo-cell used for normalization. 
# Written as a separate function so that bplapply operates in the scran namespace.
{
    if (is_sparse(x)) {
        x <- as(x, "dgCMatrix")
    } else {
        x <- as.matrix(x)
    }

    if (is.null(scaling)) {
        scaling <- colSums(x)
    }
    if (any(scaling==0)) {
        stop("cells should have non-zero library sizes or 'scaling' values")
    }
    exprs <- normalizeCounts(x, size.factors=scaling, center.size.factors=FALSE, log=FALSE)

    ave.cell <- rowMeans(exprs) * mean(scaling) # equivalent to calculateAverage().
    high.ave <- min.mean <= ave.cell 
    use.ave.cell <- ave.cell
    if (!all(high.ave)) { 
        exprs <- exprs[high.ave,,drop=FALSE]
        use.ave.cell <- use.ave.cell[high.ave]
    }

    # Using our summation approach.
    sphere <- .generateSphere(scaling)
    sizes <- sizes[sizes <= ncol(exprs)]
    new.sys <- .create_linear_system(exprs, use.ave.cell, sphere, sizes) 
    design <- new.sys$design
    output <- new.sys$output

    # Weighted least-squares.
    QR <- qr(design)
    final.nf <- qr.coef(QR, output)
    final.nf <- final.nf * scaling

    if (any(final.nf <= 0)) {
        warning("encountered non-positive size factor estimates")
        if (positive) {
            num.detected <- colSums(exprs > 0)
            final.nf <- cleanSizeFactors(final.nf, num.detected) 
        }
    }

    list(final.nf=final.nf, ave.cell=ave.cell)
}

.generateSphere <- function(lib.sizes) 
# Sorts cells by their library sizes, and generates an ordering vector
# to arrange cells in a circle based on increasing/decreasing lib size.
{
    nlibs <- length(lib.sizes)
    o <- order(lib.sizes)
    even <- seq(2,nlibs,2)
    odd <- seq(1,nlibs,2)
    out <- c(o[odd], rev(o[even]))
    c(out, out)
}

LOWWEIGHT <- 0.000001

#' @importFrom Matrix sparseMatrix
.create_linear_system <- function(cur.exprs, ave.cell, sphere, pool.sizes) 
# Does the heavy lifting of computing pool-based size factors 
# and creating the linear system out of the equations for each pool.
{
    row.dex <- col.dex <- output <- vector("list", 2L)

    # Creating the linear system with the requested pool sizes.
    out <- pool_size_factors(cur.exprs, ave.cell, sphere - 1L, pool.sizes)
    row.dex[[1]] <- out[[1]] + 1L
    col.dex[[1]] <- out[[2]] + 1L
    output[[1]]<- out[[3]]

    # Adding extra equations to guarantee solvability.
    cur.cells <- ncol(cur.exprs)
    row.dex[[2]] <- seq_len(cur.cells) + cur.cells * length(pool.sizes)
    col.dex[[2]] <- seq_len(cur.cells)
    output[[2]] <- rep(sqrt(LOWWEIGHT) / sum(ave.cell), cur.cells) # equivalent to library size factors for each cell, but downweighted.

    # Setting up the entries of the LHS matrix.
    eqn.values <- rep(c(1, sqrt(LOWWEIGHT)), lengths(row.dex))

    # Constructing a sparse matrix.
    row.dex <- unlist(row.dex)
    col.dex <- unlist(col.dex)
    output <- unlist(output)
    design <- sparseMatrix(i=row.dex, j=col.dex, x=eqn.values, dims=c(length(output), cur.cells))

    return(list(design=design, output=output))
}

#' @importFrom stats median
#' @importFrom S4Vectors wmsg
.rescale_clusters <- function(mean.prof, ref.col, min.mean) 
# Chooses a cluster as a reference and rescales all other clusters to the reference,
# based on the 'normalization factors' computed between pseudo-cells.
{
    if (is.character(ref.col)) {
        ref.col <- which(names(mean.prof)==ref.col)
        if (length(ref.col)==0L) { 
            stop("'ref.clust' not in 'clusters'")
        }
    }

    nclusters <- length(mean.prof)
    rescaling <- numeric(nclusters)
    for (clust in seq_len(nclusters)) { 
        ref.prof <- mean.prof[[ref.col]]
        cur.prof <- mean.prof[[clust]] 

        # Filtering based on the mean of the per-cluster means (requires scaling for the library size).
        # Effectively equivalent to 'calculateAverage(cbind(ref.ave.count, cur.ave.count))' where the averages
        # are themselves equivalent to 'calculateAverage()' across all cells in each cluster.
        cur.libsize <- sum(cur.prof)
        ref.libsize <- sum(ref.prof)
        to.use <- (cur.prof/cur.libsize + ref.prof/ref.libsize)/2 * (cur.libsize + ref.libsize)/2 >= min.mean
        if (!all(to.use)) { 
            cur.prof <- cur.prof[to.use]
            ref.prof <- ref.prof[to.use]
        } 

        # Adjusting for systematic differences between clusters.
        rescale.sf <- median(cur.prof/ref.prof, na.rm=TRUE)
        if (!is.finite(rescale.sf) || rescale.sf <= 0) {
            warning(wmsg("inter-cluster rescaling factor for cluster ", clust, 
                " is not strictly positive, reverting to the ratio of average library sizes"))
            rescale.sf <- sum(cur.prof)/sum(ref.prof)
        }

        rescaling[[clust]] <- rescale.sf
    }

    names(rescaling) <- names(mean.prof)
    rescaling
}

.limit_cluster_size <- function(clusters, max.size) 
# Limits the maximum cluster size to avoid problems with memory in Matrix::qr().
# Done by arbitrarily splitting large clusters so that they fall below max.size.
{
    if (!is.null(max.size) && any(table(clusters) > max.size)) { 
        clusters <- as.character(clusters)

        # NOTE: we must append '-1', even to the clusters that fall below the
        # max.size, so as to avoid name conflicts, e.g., if one cluster was
        # called "A-1" and another was called "A", appending "-1" to the latter
        # but not the former would cause issues.
        for (id in unique(clusters)) {
            current <- id==clusters
            ncells <- sum(current)
            mult <- ceiling(ncells/max.size)
            realloc <- rep(seq_len(mult), length.out=ncells)
            clusters[current] <- sprintf("%s-%s", id, realloc)
        }
    }

    clusters
}

#############################################################
# S4 method definitions.
#############################################################

#' @export
#' @rdname computePooledFactors
setGeneric("pooledSizeFactors", function(x, ...) standardGeneric("pooledSizeFactors"))

#' @export
#' @rdname computePooledFactors
setMethod("pooledSizeFactors", "ANY", .calculate_pooled_factors)

#' @export
#' @rdname computePooledFactors
#' @importFrom SummarizedExperiment assay
setMethod("pooledSizeFactors", "SummarizedExperiment", function(x, ..., assay.type="counts") {
    .calculate_pooled_factors(assay(x, i=assay.type), ...)
})

#' @export
#' @rdname computePooledFactors
#' @importFrom SummarizedExperiment assay 
#' @importFrom BiocGenerics "sizeFactors<-"
computePooledFactors <- function(x, ..., assay.type="counts") {
    sizeFactors(x) <- .calculate_pooled_factors(assay(x, i=assay.type), ...) 
    x
}
