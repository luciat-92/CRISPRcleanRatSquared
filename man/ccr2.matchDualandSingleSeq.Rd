% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISPRcleanRatSquared.R
\name{ccr2.matchDualandSingleSeq}
\alias{ccr2.matchDualandSingleSeq}
\title{Match sequences between sgRNA in dual and single libraries}
\usage{
ccr2.matchDualandSingleSeq(dual_library, single_library)
}
\arguments{
\item{dual_library}{Data frame containing information from the dual library (one row per guide pair).}

\item{single_library}{Data frame containing information from the single library (one row per guide).}
}
\value{
A data frame with one row per sgRNA in dual library (union of position 1 and position 2), including matched single library IDs and sequences.
}
\description{
This function matches sequences between a dual library and a single library in a CRISPRcleanR^2 analysis.
It creates a data frame with all sgRNAs from dual library (no repetitions, union of those in position 1 and position2) and
includes the matched single library IDs and sequences, if a match exists.
}
