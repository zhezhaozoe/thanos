% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{barplot_depths_by_sample}
\alias{barplot_depths_by_sample}
\title{Plots aggregated abundance by sample in barplot}
\usage{
barplot_depths_by_sample(ps, fill = NULL, position = "stack", ...)
}
\arguments{
\item{ps}{A phyloseq object or a list of such objects. If a list is
provided, it merges the data from all objects.}

\item{fill}{The variable by which to fill the bars. If \code{NULL}, bars
won't be filled according to a variable. The default is \code{NULL}.}

\item{position}{The position adjustment of bars. Can be "stack",
"fill", or others as defined by \code{ggplot2}. Default is "stack".}

\item{...}{Additional arguments passed to \code{geom_col}.}
}
\value{
A \code{ggplot} object representing the bar plot.
}
\description{
This function plots the aggregated abundance of genes by sample and
allows for customization of the filling color of the bars. Unlike
\code{barplot_depths()}, this function retains the sample data, enabling
the use of \code{facet_wrap(~whatever)} for further customization.
}
\examples{
## Assuming 'ps' is a phyloseq object and 'ps_list' is a list of phyloseq objects
barplot_depths_by_sample(ps, fill = "Kingdom")
barplot_depths_by_sample(ps_list, fill = "Phylum", position = "fill")

}
