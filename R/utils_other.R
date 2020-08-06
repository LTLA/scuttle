.replace <- function(new, old) {
    if (!is.null(old)) old else new
}

.use_names_to_integer_indices <- function(use, x, nameFUN, msg) {
    if (isTRUE(use)) {
        use <- seq_along(nameFUN(x))
    } else if (isFALSE(use)) {
        use <- integer(0)
    } else if (is.character(use)) {
        use <- match(use, nameFUN(x))
        if (any(is.na(use))) {
            stop(sprintf("'%s' contains invalid values", msg))
        }
    } else {
        if (any(use < 1L) || any(use > length(nameFUN(x)))) {
            stop(sprintf("'%s' contains out-of-bounds indices", msg))
        }
    }
    use
}

#' @export
#' @importClassesFrom S4Vectors List DataFrame
.unpackLists <- function(...) {
    objects <- list(...)
    for (i in seq_along(objects)) {
        current <- objects[[i]]
        if (!is.list(current)) {
            if (is(current, "List") && !is(current, "DataFrame")) {
                current <- as.list(current)
            } else {
                current <- list(current)
            }
            objects[[i]] <- current
        }
    }
    do.call(c, objects)
}
