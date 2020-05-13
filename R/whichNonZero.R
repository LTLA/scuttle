#' Find non-zero entries of a matrix
#'
#' Finds the non-zero entries of a matrix in the most efficient manner.
#' Not sure there's much more to say here.
#' 
#' @param x A numeric matrix-like object, usually sparse in content if not in representation.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object controlling how parallelization should be performed.
#' only used when \code{x} is a \linkS4class{DelayedMatrix} object.
#' @param ... For the generic, additional arguments to pass to the specific methods.
#'
#' For the methods, additional arguments that are currently ignored.
#'
#' @return A list containing \code{i}, an integer vector of the row indices of all non-zero entries;
#' \code{j}, an integer vector of the column indices of all non-zero entries;
#' and \code{x}, a numeric vector of the values of the non-zero entries.
#'
#' @author Aaron Lun
#'
#' @examples
#' x <- Matrix::rsparsematrix(1e6, 1e6, 0.000001)
#' out <- whichNonZero(x)
#' str(out)
#' 
#' @seealso
#' \code{\link{which}}, obviously.
#'
#' @export
setGeneric("whichNonZero", function(x, ...) standardGeneric("whichNonZero"))

#' @export
#' @rdname whichNonZero
#' @importClassesFrom Matrix dgTMatrix
setMethod("whichNonZero", "dgTMatrix", function(x, ...) {
    list(i=x@i+1L, j=x@j+1L, x=x@x)
})

#' @export
#' @rdname whichNonZero
#' @importClassesFrom Matrix dgCMatrix
setMethod("whichNonZero", "dgCMatrix", function(x, ...) {
    i <- x@i + 1L
    d <- diff(x@p)
    j <- rep(seq_along(d), d)
    list(i=i, j=j, x=x@x)
})

#' @export
#' @rdname whichNonZero
#' @importFrom Matrix which
setMethod("whichNonZero", "ANY", function(x, ...) {
    keep <- which(x!=0, arr.ind=TRUE)
    list(i=keep[,1], j=keep[,2], x=x[keep])
})

#' @export
#' @rdname whichNonZero
#' @importFrom DelayedArray blockApply getAutoBPPARAM setAutoBPPARAM
#' @importFrom BiocParallel SerialParam
setMethod("whichNonZero", "DelayedArray", function(x, BPPARAM=SerialParam(), ...) {
    oldBP <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldBP))

    output <- blockApply(x, FUN=.find_blocked_nzero)
    do.call(mapply, c(list(FUN=c), output, list(SIMPLIFY=FALSE)))
})

#' @importFrom Matrix which
#' @importFrom IRanges start
.find_blocked_nzero <- function(block) {
    grid <- attr(block, "from_grid")
    bid <- attr(block, "block_id") 
    viewport <- grid[[bid]]

    basic <- which(block!=0, arr.ind=TRUE) # guaranteed to be ordinary, no need for Matrix::which.
    i <- unname(basic[,1]) + as.integer(start(viewport)[1]) - 1L
    j <- unname(basic[,2]) + as.integer(start(viewport)[2]) - 1L

    list(i=i, j=j, x=block[basic])
}
