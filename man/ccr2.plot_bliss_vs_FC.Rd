% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISPRcleanRatSquared.R
\name{ccr2.plot_bliss_vs_FC}
\alias{ccr2.plot_bliss_vs_FC}
\title{Plot Bliss Z-score vs logFC to identify synergistic and lethal pairs.}
\usage{
ccr2.plot_bliss_vs_FC(
  dual_FC,
  corrected = FALSE,
  THR_FC = -1,
  THR_BLISS = -1,
  saveToFig = FALSE,
  saveFormat = "pdf",
  outdir = "./",
  EXPname = ""
)
}
\arguments{
\item{dual_FC}{A data frame containing information about dual KO logFCs.
It should have columns 'avgFC_scaled', 'correctedFC_scaled' (see ccr2.scale_pos_neg), and 'bliss_zscore',
'bliss_zscore_corrected' (see ccr2.compute_bliss), and 'info'.}

\item{corrected}{Logical, indicating whether to use CRISPRcleanR^2 corrected values for logFC and Bliss Z-score (default is FALSE).}

\item{THR_FC}{Numeric, the threshold for logFC below which points will be considered lethal and stored (default is -1).}

\item{THR_BLISS}{Numeric, the threshold for Bliss Z-score below which points will be considered synergistic and stored (default is -1).}

\item{saveToFig}{Logical, indicating whether to save the plot to a file (default is FALSE).}

\item{saveFormat}{Character, the format in which to save the plot (default is "pdf").}

\item{outdir}{Character, the directory where the plot file will be saved (default is "./").}

\item{EXPname}{Character, an additional name to include in the plot file name (default is "").}
}
\value{
A data frame containing dual guide pairs classified as lethal and synergistic based on the specified thresholds. Same structure as input dual_FC data frame.
}
\description{
This function generates a scatter plot comparing Bliss Z-score with logFC in dual KO screen.
Here, scaled logFCs are used, centering positive controls to -1 and negative controls to 0 (see ccr2.scale_pos_neg).
It can be used to plot results both before and after CRISPRcleanR^2 correction.
Dual guides with bliss z-score and logFCs passing pre-specified thresholds are stored.
}
\seealso{
\code{\link{ccr2.compute_bliss}},
\code{\link{ccr2.scale_pos_neg}}
}
