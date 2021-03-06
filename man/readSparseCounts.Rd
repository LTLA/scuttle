% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readSparseCounts.R
\name{readSparseCounts}
\alias{readSparseCounts}
\title{Read sparse count matrix from file}
\usage{
readSparseCounts(
  file,
  sep = "\\t",
  quote = NULL,
  comment.char = "",
  row.names = TRUE,
  col.names = TRUE,
  ignore.row = 0L,
  skip.row = 0L,
  ignore.col = 0L,
  skip.col = 0L,
  chunk = 1000L
)
}
\arguments{
\item{file}{A string containing a file path to a count table, or a connection object opened in read-only text mode.}

\item{sep}{A string specifying the delimiter between fields in \code{file}.}

\item{quote}{A string specifying the quote character, e.g., in column or row names.}

\item{comment.char}{A string specifying the comment character after which values are ignored.}

\item{row.names}{A logical scalar specifying whether row names are present.}

\item{col.names}{A logical scalar specifying whether column names are present.}

\item{ignore.row}{An integer scalar specifying the number of rows to ignore at the start of the file, \emph{before} the column names.}

\item{skip.row}{An integer scalar specifying the number of rows to ignore at the start of the file, \emph{after} the column names.}

\item{ignore.col}{An integer scalar specifying the number of columns to ignore at the start of the file, \emph{before} the column names.}

\item{skip.col}{An integer scalar specifying the number of columns to ignore at the start of the file, \emph{after} the column names.}

\item{chunk}{A integer scalar indicating the chunk size to use, i.e., number of rows to read at any one time.}
}
\value{
A dgCMatrix containing double-precision values (usually counts) for each row (gene) and column (cell).
}
\description{
Reads a sparse count matrix from file containing a dense tabular format.
}
\details{
This function provides a convenient method for reading dense arrays from flat files into a sparse matrix in memory.
Memory usage can be further improved by setting \code{chunk} to a smaller positive value.

The \code{ignore.*} and \code{skip.*} parameters allow irrelevant rows or columns to be skipped.
Note that the distinction between the two parameters is only relevant when \code{row.names=FALSE} (for skipping/ignoring columns) or \code{col.names=FALSE} (for rows).
}
\examples{
outfile <- tempfile()
write.table(data.frame(A=1:5, B=0, C=0:4, row.names=letters[1:5]), 
    file=outfile, col.names=NA, sep="\t", quote=FALSE)

readSparseCounts(outfile)

}
\seealso{
\code{\link{read.table}},
\code{\link{readMM}}
}
\author{
Aaron Lun
}
