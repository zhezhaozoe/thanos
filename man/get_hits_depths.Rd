% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmmer.R
\name{get_hits_depths}
\alias{get_hits_depths}
\title{Compare OTU hits depths between query and control}
\usage{
get_hits_depths(
  ps,
  query_tblout,
  control_tblout,
  linker,
  taxrank = NULL,
  phyloseq = TRUE
)
}
\arguments{
\item{ps}{A phyloseq object containing OTU (Operational Taxonomic
Units) counts and (optionally) taxonomic information.}

\item{query_tblout}{Character; the HMMER output for the query genes.}

\item{control_tblout}{Character; the HMMER output for the control
genes.}

\item{linker}{Function that links OTU identifiers between tblout
tables and the phyloseq object. The function must take two arguments:
SeqFile and Target.}

\item{taxrank}{Optional; character specifying the taxonomic rank
at which to aggregate hits. Valid options depend on the taxonomic
ranks present in the \code{ps} object. If not specified, no aggregation is
performed.}

\item{phyloseq}{Logical; if TRUE, returns a phyloseq object,
otherwise returns a data.table object. Defaults to TRUE.}
}
\value{
If \code{phyloseq} is TRUE (the default), this function returns
a phyloseq object that contains the relative OTU depths of query vs
control genes (aggregated by taxrank if specified). If \code{phyloseq} is
FALSE, it returns a data.table object containing the calculated depths
for queries and controls separately.
}
\description{
\code{get_hits_depths} function calculates the relative OTU
depths between the query and control genes, optionally aggregating hits
at a specified taxonomic rank.
}
\details{
The function prunes the phyloseq object to include only taxa
that are present in both the query and control datasets based on the
\code{linker} function. If a taxonomic rank is specified, the function
aggregates hits at that rank and recalculates taxa names. It also
warns about queries not found in the control. It performs calculates
the depth ratios between query and control conditions. Two linker
functions are provided: mags_linker() and contigs_linker().
}
