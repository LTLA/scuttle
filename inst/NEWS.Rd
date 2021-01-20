\name{NEWS}
\title{News for Package \pkg{scuttle}}

\section{Version 1.2.0}{\itemize{
\item Migrated \code{whichNonZero()} to \pkg{beachmat}.

\item Bugfixes for factor-based \code{colData} aggregation in \code{aggregateAcrossCells()}.
Added proper support for Vectors.

\item Bugfix for correct response to \code{use.altexps=} in \code{perCellQCMetrics()}, \code{perFeatureQCMetrics()}.

\item Added a \code{normalize.all=} option to \code{normalizeCounts()}.
Removed unnecessary warning when \code{down.target=} is not specified.
Exposed the default \code{size.factors=} in the SingleCellExperiment method.

\item Modified the SingleCellExperiment method of \code{logNormCounts()} so that manually specified size factors do not apply to alternative Experiments.
Only relevant if \code{size.factors=} and \code{use.altexps=} are specified.
}}

\section{Version 1.0.0}{\itemize{
\item Split off \pkg{scuttle} from \pkg{scater} by migrating all non-visualization code from the latter.

\item Began transition to dot-separate argument names from original snake case format.

\item Added a \code{geometricSizeFactors()} function, deprecated \code{geometric=TRUE} in \code{librarySizeFactors()}.

\item Single-object downsampling in \code{downsampleBatches()} now behaves more consistently with multi-object downsampling.
}}

