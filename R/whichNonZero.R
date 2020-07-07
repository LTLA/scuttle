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
    j <- rep.int(seq_along(d), d)
    list(i=i, j=j, x=x@x)
})

#' @export
#' @rdname whichNonZero
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM nzindex nzdata
#' @importClassesFrom DelayedArray SparseArraySeed
#' @importFrom BiocParallel SerialParam
setMethod("whichNonZero", "ANY", function(x, BPPARAM=SerialParam(), ...) {
    oldBP <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldBP))

    out <- as(x, "SparseArraySeed")
    idx <- nzindex(out)
    list(i=idx[,1], j=idx[,2], x=nzdata(out))
})
