With the ever growing amount of metagenomic DNA sequences, the need of efficient and robust methods for functional analysis of the samples also grows.
Here, we present thanos, an R package for the analysis and visualisation of gene functions in metagenomic samples.
The package aims to be a low-code solution to a common problem in environmental metagenomics: finding the prevalence of pathways in a sample, using the count of genes as a proxy.

The central entity in thanos is the depths matrix, which contains the sequencing depths of the elements of interest.
The elements can be either assembled contigs, binned metagenomes, or general operational taxonomic units (OTUs), including marker OTUs (mOTUs, @ref).

# importing files
The package leverages and extends the `phyloseq` package, which already provides a good foundation for common operations on depths files.
Utility functions are provided to import each type of file directly from its path, returning an otu_table() with a consistent interface.
In the case of contigs, the raw depths are normalised as $contig_depth * contig_length / sum(total_contig_depth * contig_length)$.
There are also functions to import taxonomy tables from the widely used  GTDB toolkit results.

# searching genes of interest
The main function of the package, though, is to perform gene searches and analyse the results.
Users can start from the Kegg Ortholog number of their gene of interest (say, K432148 for hzs1), get an HMM profile for that gene based on all the kegg orthologs, and search the profile in the gene sequences of the elements of interest (either contigs or MAGs).
By matching the search hits with the elements' rownames, we are able to calculate the total depth of the gene of interest.
However, this alone would not be enough without a control.
Thanos uses universally conserved single-copy genes as a normalisation control for the depth of the genes of interest.
Briefly, the depth of the gene of interest is divided by the depth of the single-copy genes, providing an estimate of the relative copy number of the gene of interest.
The analysis can be performed either in bulk for the whole samples, or conditionally by taxa.

In the case of contig depths, whose files typically have a considerable size, we provide an alternative workflow that is memory-efficient.

# visualising the results
Good visualisations are a pivotal tool for data analysis.
For this reason, thanos exposes functions to plot functional or compositional profiles.
Thanks to the pow

# customising the plots
An important part of achiving beautiful visualization is the ability to choose the order of the categorical variables on the x axis.
The package provides two filter function `sort_samples_by()`, and `set_group_order()`, which return a modified phyloseq object which stores information about the ordering of samples or other categorical values.
Similarly, users may want to customise the colors in the plot.
This can be accomplished with the `...` function.
