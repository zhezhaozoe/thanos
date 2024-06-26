% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{read_gtdbtk}
\alias{read_gtdbtk}
\title{Read GTDB-Tk classification output}
\usage{
read_gtdbtk(file)
}
\arguments{
\item{file}{A string specifying the path to the GTDB-Tk
classification output file.}
}
\value{
A matrix where rows are genomes (user_genome) and columns are
the taxonomic categorizations (Domain, Phylum, Class, Order, Family,
Genus, Species).
}
\description{
Reads a GTDB-Tk (Genome Taxonomy Database Toolkit) classification
file and extracts the taxonomic classification into a
phyloseq-friendly format. It parses the classifications into domain,
phylum, class, order, family, genus, and species, and returns a
matrix with these categories as column names and genomes as row
names.
}
\examples{
# Assume "gtdbtk_output.tsv" is a file in the working directory with GTDB-Tk classifications
gtdbtk_results <- read_gtdbtk("gtdbtk_output.tsv")
head(gtdbtk_results)

}
