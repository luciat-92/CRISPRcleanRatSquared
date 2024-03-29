% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISPRcleanRatSquared.R
\name{ccr2.run_complete}
\alias{ccr2.run_complete}
\title{Run complete CRISPRcleanR^2 analysis for dual screens.}
\usage{
ccr2.run_complete(
  filename_single,
  min_reads_single = 30,
  libraryAnnotation_single,
  min_reads = 30,
  libraryAnnotation_dual,
  dual_count = NULL,
  dual_logFC = NULL,
  EXPname,
  display = FALSE,
  outdir = "./",
  saveToFig,
  saveFormat = "pdf",
  correctGW,
  excludeGene_plot = NULL,
  CNA,
  CN_thr = 8,
  ...
)
}
\arguments{
\item{filename_single}{A string specifying the path of a tsv file containing the raw sgRNA counts. This must be a tab delimited file with one row per sgRNA and the following columns/headers:
\itemize{
\item sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
\item gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration;
followed by the columns containing the sgRNAs' counts for the controls and columns for library trasfected samples.
}}

\item{min_reads_single}{This parameter defines a filter threshold value for sgRNAs, based on their average counts in the control sample.
Specifically, it indicates the minimal number of counts that each individual sgRNA needs to have in the controls (on average) in order to be included in the output.}

\item{libraryAnnotation_single}{A data frame containing the sgRNA annotations, with a named row for each sgRNA, and columns for targeted genes, genomic coordinates and possibly other information.}

\item{min_reads}{Minimum number of reads required in the control sample (e.g. plasmid) for a pair to be retained (default is 30).}

\item{libraryAnnotation_dual}{Data frame containing information from the dual library (one row per guide pair).}

\item{dual_count}{Data frame containing the raw counts from dual KO screen, raw counts start from the 9th column. The 9th column include the control raw counts (default is NULL).}

\item{dual_logFC}{Data frame containing dual KO log fold changes (default is NULL).}

\item{EXPname}{Name of the experiment.}

\item{display}{Logical, whether to display plots (default is FALSE).}

\item{outdir}{Directory to save the output files (default is "./").}

\item{saveToFig}{Logical, whether to save figures.}

\item{saveFormat}{File format for saving figures (default is "pdf").}

\item{correctGW}{A character specifying the type of pseudo-single logFCs mean centering after data-injection ("CHR": per chromosome, "GW": genome-wide, or NULL for no centering). Default is NULL. #TODO: USED FOR TESTING, REMOVE!}

\item{CNA}{A data frame containing copy number alteration information. It must include the following columns:
\itemize{
\item CHROM: chromosome id in "chri" format
\item start: start of the CN region, genomic location in base pair
\item end: end of the CN region, genomic location in base pair
\item C: copy number (absolute) value
}}

\item{CN_thr}{Numeric (absolute), a threshold for defining extreme CNA.}

\item{...}{Additional parameters to be passed to other functions.}

\item{excludeGene}{Character vector, genes to be excluded from dual_FC_correctedFC (default is NULL).}
}
\value{
A list containing:
\itemize{
\item dual: data frame of dual KO logFCs including CRISPRcleanR^2 corrected values, bliss z-score synergy, scaled logFCs per positive and negative controls (one row per guide pair),
\item dual_gene: data frame same as dual but at the gene pair level (scores collapsed using median),
\item single_gw: data frame of single KO logFCs including CRISPRcleanR corrected values at the genome-wide level,
\item single: data frame of single KO logFCs including CRISPRcleanR corrected values only for those sgRNAs also present in dual KO library,
\item pseudo_single_segments: identified genomic regions (segments) by CRISPRcleanR applied to pseudo single logFCs,
\item system_solition: data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides (both positions), it also includes the fitted pseudo-single correction obtained as system solution,
\item model_perf:  model performance for single VS pseudo-single, including both double-cut pairs and non-target pairs,
\item model_est: model estimates (coefficients) for single VS pseudo-single, including both double-cut pairs and non-target pairs,
\item top_corrected: using CRISPRcleanR^2 corrected results, dual guide pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
\item top_gene_corrected: using CRISPRcleanR^2 corrected results, dual gene pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
\item top_uncorrected: using uncorrected (original) results, dual guide pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
\item top_gene_uncorrected: using uncorrected (original) results, dual gene pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
\item max_correction: maximum dual correction allowed based on single KO correction (used as upper bound for the optimization problem solution).,
\item min_correction: minimum dual correction allowed based on single KO correction (used as upper bound for the optimization problem solution).
}
}
\description{
This function is the complete wrap-up for the entire pipeline, from raw count of both single and dual KO screens to CRISPRcleanR^2 correction, synergy estimation and results visualization.
In particular, this function
\itemize{
\item runs CRISPRcleanR pipeline on single KO raw counts (genome-wide) (see \code{\link{ccr.run_complete}}),
\item matches single KO and dual KO libraries based on sequence (see \code{\link{ccr2.matchDualandSingleSeq}}),
\item pre-processes dual KO screen, converting raw counts into logFCs (if necessary, see \code{\link{ccr2.NormfoldChanges}}), averages logFCs across replicates and annotate sgRNA pairs with chromosomal positions (see \code{\link{ccr2.logFCs2chromPos}}),
\item computes singletons (guide pairs matched with non-essential gene or intergenic region) summary statistics to be used in bliss model (see \code{\link{ccr2.get_summary_singletons}}),
\item computes dual KO logFCs correction for "non-target" pairs (see \code{\link{ccr2.run_nontarget}}),
\item computes dual KO logFCs correction for "double-cuts" dual pairs, targeting two existing genomic regions (see \code{\link{ccr2.run}}),
\item combines the results from both non-target and double-cuts pairs,
\item computes synergy based on the bliss model for both uncorrected and corrected dual KO logFCs (see \code{\link{ccr2.compute_bliss}}),
\item visualize results such as relationship with CNA status (see \code{\link{ccr2.plot_correction}}, \code{\link{ccr2.plotClasses}}, \code{\link{ccr2.plotCNA}}, \code{\link{ccr2.plotCNAdensity}}),
\item visualize results for single KO for guides available also in with dual KO library (see \code{\link{ccr2.plotMatchingSingle}}),
\item centers logFCs of positive controls to -1 and negative controls to 0 for both uncorrected and corrected logFCs (see \code{\link{ccr2.scale_pos_neg}}),
\item selects synergistic and lethal guide pairs based on logFCs threshold -1 and bliss z-score threshold 0.5 (see \code{\link{ccr2.plot_bliss_vs_FC}}),
\item compute logFCs at the gene pair level (median) and repeats the last 2 steps to obtain synergistic and lethal gene pairs.
}
}
\seealso{
\code{\link{ccr.run_complete}},
\code{\link{ccr2.matchDualandSingleSeq}},
\code{\link{ccr2.NormfoldChanges}},
\code{\link{ccr2.logFCs2chromPos}},
\code{\link{ccr2.get_summary_singletons}},
\code{\link{ccr2.match_singletons_singleFC}},
\code{\link{ccr2.run_nontarget}},
\code{\link{ccr2.run}},
\code{\link{ccr2.compute_bliss}},
\code{\link{ccr2.plot_correction}},
\code{\link{ccr2.plotClasses}},
\code{\link{ccr2.plotCNA}},
\code{\link{ccr2.plotCNAdensity}},
\code{\link{ccr2.plotMatchingSingle}},
\code{\link{ccr2.scale_pos_neg}},
\code{\link{ccr2.plot_bliss_vs_FC}}
}
