% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISPRcleanRatSquared.R
\name{ccr2.filterGWclean}
\alias{ccr2.filterGWclean}
\title{Post-processing of CRISPRcleanR corrected psuedo-single logFCs}
\usage{
ccr2.filterGWclean(
  dataInjection_correctedFCs,
  dataInjection_segments,
  pseudo_single,
  guide_id,
  saveToFig = FALSE,
  display = TRUE,
  saveFormat = ".pdf",
  outdir = "./",
  EXPname = ""
)
}
\arguments{
\item{dataInjection_correctedFCs}{Data frame containing corrected injected pseudo-single logFCs via CRISPRcleanR (output of ccr.GWclean).}

\item{dataInjection_segments}{Data frame containing information about identified gene segments from injected pseudo-single logFCs via CRISPRcleanR (output of ccr.GWclean).}

\item{pseudo_single}{Data frame containing information for pseudo-single logFCs (before CRISPRcleanR correction, output of ccr2.createPseudoSingle_combine)}

\item{guide_id}{Identifier for the guide position (1 or 2).}

\item{saveToFig}{Logical, indicating whether to save the generated plots to files. Default is FALSE.}

\item{display}{Logical, indicating whether to display the plots. Default is TRUE.}

\item{saveFormat}{File format for saving the plots (e.g., ".pdf", ".png"). Default is ".pdf".}

\item{outdir}{Directory to save the plots. Default is "./".}

\item{EXPname}{Experiment name to include in the plot filenames.}
}
\value{
A list containing data frames with information about corrected pseudo-single logFCs
and gene segments identified by CRISPRcleanR, for those guides originally available in the dual KO screen.
}
\description{
This function filters and post-process psuedo-single logFCs after genome-wide
CRISPRcleanR correction (output of ccr.GWclean). It filters for those guides originally available in the pseudo-single logFCs (i.e. present in dual KO)
It provides information about the pseudo-single correction values and gene segments in a specific guide position (1 or 2).
}
