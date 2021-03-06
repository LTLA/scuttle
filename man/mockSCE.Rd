% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mockSCE.R
\name{mockSCE}
\alias{mockSCE}
\title{Mock up a SingleCellExperiment}
\usage{
mockSCE(ncells = 200, ngenes = 2000, nspikes = 100)
}
\arguments{
\item{ncells}{Integer scalar, number of cells to simulate.}

\item{ngenes}{Integer scalar, number of genes to simulate.}

\item{nspikes}{Integer scalar, number of spike-in transcripts to simulate.}
}
\value{
A SingleCellExperiment object containing a count matrix in the \code{"counts"} assay,
a set of simulated \code{\link{colData}} fields,
and spike-in data in the \code{"Spikes"} field of \code{\link{altExps}}.
}
\description{
Mock up a \linkS4class{SingleCellExperiment} containing simulated data,
for use in documentation examples.
}
\details{
Users should set a seed to obtain reproducible results from this function.
}
\examples{
set.seed(1000)
sce <- mockSCE()
sce

}
\seealso{
\code{\link{SingleCellExperiment}}, for the constructor.
}
\author{
Aaron Lun
}
