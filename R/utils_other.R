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
