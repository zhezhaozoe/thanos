% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shortcuts.R
\name{build_hmm_from_ko}
\alias{build_hmm_from_ko}
\title{Build an HMM profile for the given KO identifiers}
\usage{
build_hmm_from_ko(kos, method = "Muscle", hmmer_path = "", ...)
}
\arguments{
\item{method}{Character. The method used for performing the MSA.
Defaults to "Muscle".}

\item{...}{Additional arguments passed to kegg_kegg_msa().}

\item{ko}{A character string representing the KEGG Orthology (KO)
identifier.}
}
\value{
The path to the HMM file.
}
\description{
Given a KO identifier, this function retrieves the corresponding
protein sequences and builds a Multiple Sequence Alignment (MSA), then
builds a Hidden Markov Model (HMM) based on the MSA.
}
\examples{
hmm_from_ko("K00001")

}
