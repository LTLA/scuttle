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

\item Deprecated \code{use.altexps=} in favor of \code{applySCE()} in \code{logNormCounts()} and \code{aggregateAcrossCells()}.

\item Renamed \code{addPerCellQC()} and \code{addPerFeatureQC()} to \code{addPerCellQCMetrics()} and \code{addPerCellFeatureMetrics()}, for consistency.
Soft-deprecated the old functions.

\item Moved most of \code{quickPerCellQC()} functionality into the new \code{perCellQCFilters()} function.
Repurposed the former to directly return a filtered SummarizedExperiment object.

\item Migrated \pkg{scran}'s normalization-related functions into this package.
Added \code{pooledSizeFactors()}, \code{computePooledFactors()}, \code{cleanSizeFactors()} and \code{computeSpikeFactors()}.

\item Added \code{transform="asinh"} to \code{normalizeCounts()} and \code{logNormCounts()} for inverse hyperbolic transformations of CITE-seq data.
}}

\section{Version 1.0.0}{\itemize{
\item Split off \pkg{scuttle} from \pkg{scater} by migrating all non-visualization code from the latter.

\item Began transition to dot-separate argument names from original snake case format.

\item Added a \code{geometricSizeFactors()} function, deprecated \code{geometric=TRUE} in \code{librarySizeFactors()}.

\item Single-object downsampling in \code{downsampleBatches()} now behaves more consistently with multi-object downsampling.
}}

