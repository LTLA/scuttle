#' Single-cell utilities 
#'
#' @description
#' The \pkg{scuttle} package provides some utility functions for single-cell 'omics data analysis.
#' This includes some simple methods for computing and filtering on quality control;
#' basic data transformations involving various types of scaling normalization;
#' and flexible aggregation across cells or features.
#'
#' \pkg{scuttle} also implements wrapper functions that simplify boilerplate for developers of client packages.
#' This includes packages such as \pkg{scran}, \pkg{scater} and \pkg{DropletUtils}, to name a few.
#' Note that much of the code here was inherited from the \pkg{scater} package.
#'
#' @author Aaron Lun
#' 
#' @name scuttle-pkg
#' @useDynLib scuttle, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import methods
NULL

