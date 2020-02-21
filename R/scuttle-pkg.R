#' The \pkg{scuttle} package
#'
#' The \pkg{scuttle} package provides some utility functions from single-cell 'omics data analysis.
#' It focuses on data transformations such as quality control and normalization,
#' which are prerequisites for the majority of downstream analysis pipelines.
#' It provides methods for flexible aggregation across cells or features,
#' which are typically used to obtain summary statistics for groups of cells or for gene sets.
#'
#' \pkg{scuttle} also implements wrapper functions that simplify boilerplate for developers of client packages.
#' This includes packages such as \pkg{scran}, \pkg{scater} and \pkg{DropletUtils}, to name a few.
#' Note that much of the code here was inherited from the \pkg{scater} package.
#'
#' @author Aaron Lun
#' 
#' @rdname scuttle-pkg
#' @useDynLib scuttle, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import methods
NULL

