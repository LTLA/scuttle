% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePerFeatureDF.R
\name{makePerFeatureDF}
\alias{makePerFeatureDF}
\title{Create a per-feature data.frame}
\usage{
makePerFeatureDF(
  x,
  cells = NULL,
  assay.type = "logcounts",
  use.rowdata = TRUE,
  check.names = FALSE,
  exprs_values = NULL,
  check_names = NULL
)
}
\arguments{
\item{x}{A \linkS4class{SingleCellExperiment} object.
This is expected to have non-\code{NULL} row names.}

\item{cells}{Character vector specifying the features for which to extract expression profiles across cells.}

\item{assay.type}{String or integer scalar indicating the assay to use to obtain expression values.
Must refer to a matrix-like object with integer or numeric values.}

\item{use.rowdata}{Logical scalar indicating whether row metadata of \code{x} should be included.
Alternatively, a character or integer vector specifying the row metadata fields to use.}

\item{check.names}{Logical scalar indicating whether column names of the output should be made syntactically valid and unique.}

\item{exprs_values, check_names}{Soft-deprecated equivalents to the arguments above.}
}
\value{
A data.frame containing one field per aspect of data in \code{x} - see Details.
Each row corresponds to a feature (i.e., row) of \code{x}.
}
\description{
Create a per-feature data.frame (i.e., where each row represents a feature) from a \linkS4class{SingleCellExperiment},
most typically for creating custom \pkg{ggplot2} plots.
}
\details{
This function enables us to conveniently create a per-feature data.frame from a \linkS4class{SingleCellExperiment}.
Each row of the returned data.frame corresponds to a row in \code{x},
while each column of the data.frame corresponds to one aspect of the (meta)data in \code{x}.

Columns are provided in the following order:
\enumerate{
\item Columns named according to the entries of \code{cells} represent the expression values across features for the specified cell in the \code{assay.type} assay.
\item Columns named according to the columns of \code{rowData(x)} represent the row metadata variables.
This consists of all variables if \code{use.rowdata=TRUE}, no variables if \code{use.rowdata=FALSE},
and only the specified variables if \code{use.rowdata} is set to an integer or character vector.
}

By default, nothing is done to resolve syntactically invalid or duplicated column names.
\code{check_names=TRUE}, this is resolved by passing the column names through \code{\link{make.names}}.
Of course, as a result, some columns may not have the same names as the original fields in \code{x}.
}
\examples{
example_sce <- mockSCE()
example_sce <- logNormCounts(example_sce)
rowData(example_sce)$Length <- runif(nrow(example_sce))

df <- makePerFeatureDF(example_sce, cells="Cell_001")
head(df)

}
\seealso{
\code{\link{makePerCellDF}}, for the cell-level equivalent.
}
\author{
Aaron Lun
}
