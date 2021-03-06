% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{EpiDatatype-class}
\alias{EpiDatatype-class}
\title{Class for an epigenomics data type to be added to the analysis.}
\description{
Class for an epigenomics data type to be added to the analysis.
}
\section{Slots}{

\describe{
\item{\code{name}}{The name of the epigenomic data type.}

\item{\code{table}}{A data frame with each row corresponding to a sample. Requires `experiment` and `trt.file` columns.
If \code{has.control} is \code{TURE}, `ctrl.file` column is required..}

\item{\code{experiments}}{Names of the experiments used for the analysis.}

\item{\code{file.type}}{The type of the alignment files, currently supporting .bam and .tagAlign.}

\item{\code{has.control}}{A logical value indicating if the data type has input samples.}

\item{\code{data.dir}}{The directory of data files}

\item{\code{sig.range}}{A \code{\link{GenomicRanges}} object prespecifing where to look for signals in the genome.}

\item{\code{bin.width}}{The length of bins to be used for summarizing counts.}

\item{\code{chrlen.file}}{A text file with the length of chromosomes of the target genome.}

\item{\code{chr.count}}{The number of chromosomes there is.}
}}

