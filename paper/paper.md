---
title: 'Thanos: an R pacakge for gene-centric analysis of the functional potential in metagenomic samples'
tags:
  - R
  - metagenomics
  - functional profiling
authors:
  - name: Zhe Zhao
    orcid: 0000-0001-9515-9919
    affiliation: "1, 2"
  - name: Federico Marotta
    orcid: 0000-0002-0174-3901
    affiliation: 2
  - name: Min Wu
    affiliation: 1
affiliations:
 - name: College of Life Sciences, Zhejiang University, People’s Republic of China
   index: 1
 - name: European Molecular Biology Laboratory, Germany
   index: 2

date: 25 March 2024
bibliography: paper.bib
---

# Abstract

As the amount of metagenomic sequencing keeps increasing, there is a growing need for tools that help biologists make sense of the data.
Specifically, researchers are often interested in the potential of a microbial community to carry out a metabolic reaction, but this analysis requires knitting together multiple software tools into a complex pipeline.
`thanos` offers a user-friendly R package designed for the pathway-centric analysis and visualization of the functions encoded within metagenomic samples.
It allows researchers to go beyond taxonomic profiles and find out, in a quantitative way, which pathways are prevalent in an environment, as well as comparing different environments in terms of their functional potential.
The analysis is based on the sequencing depth of the genes of interest, either in the metagenome-assembled genomes (MAGs) or in the raw data (contigs), using a normalization strategy that enables comparison across samples.
The package can import the data from multiple formats and offers functions for the visualization of the results as bar plots of the functional profile, boxplots to compare functions across samples, and annotated pathway graphs.
By streamlining the analysis of the functonal potential encoded in microbial communities, `thanos` can enable impactful discoveries in all the fields touched by metagenomics, from human health to the environmental sciences.

# Introduction

The field of metagenomics has experienced significant growth over recent decades, offering unprecedented insights into microbial communities across various environments, ranging from the human gut to marine ecosystems.
The progress was largely enabled by the advent of next-generation sequencing technologies and the development of algorithms for the reconstruction of individual microbial genomes from the pooled and fragmented DNA extracted from environmental samples.
As is often the case, however, collecting more data does not necessarily lead to an improved understaing of the underlying biological systems.
Thus, with increasing sequencing data, there also increases the need for tools capable of interpreting them.
As R is one of the most popular languages for bioinformatics, not least thanks to the Bioconductor project, and since Thanos is implemented in R, in this work we shall focus mainly on the R ecosystem.

One of the most basic questions that can be asked concerns which taxa are present in a sample.
The `phyloseq` [@mcmurdie2013phyloseq] package is a powerful tool for the exploration of microbiome profiles, offering a practical approach to investigate the taxonomic composition of metagenomic samples collected from diverse environments.
However, the taxonomic profile offers only one view into the complex and multifaceded nature of biological samples.
The gene composition of a sample offers a complementary view, one that can help answer questions such as: does carbon fixation occur in this environment?, or: is methane metabolism more active in lakes or in the Atlantic ocean?
This sort of information is not always reflected in the taxonomic profiles.
Furthermore, taxonomic profiling requires either reconstructing the full metagenome-assembled genomes (MAGs), or at least recovering several marker genes from the same operational taxonomic unit (OTU) [@motus].
On the other hand, the functional profile can be obtained on an individual gene basis, and can therefore be extracted even from unbinned contigs.
This allows researchers to overcome the "binning bias" and consider all of the reads in a sample, not just those that were binned into MAGs.

As the functional profiling task is somewhat more complex than taxonomic profiling, there is currently no standard way to perform this analysis.
Tools like Prokka [@seeman2014prokka] or the EggNog-mapper [@cantalapiedra2021eggnog] can indeed provide a bulk-level overview of the functional composition by performing sequence-similarity searches for each gene in the sample, but researchers are often interested in more specific questions about individual genes or metabolic pathways, and require a more in-depth analysis.
In this work, we present Thanos, an R package that offers a convenient way to perform functional profiling with a gene- or pathway- centric approach.
The key innovation of the package is its ability to provide quantitative functional information "depth score" for each gene of interest across samples.
Genes with a higher sequencing depth are thought to be, at least as a first approximation, more prevalent and active in the sample.
Moreover, we introduce a normalization strategy that makes depth scores comparable across samples and even across independent sequencing projects, enabling the comparative analysis of multiple environments at the same time.
The depth scores of individual genes can also be aggregated into their natural higher-level units, the metabolic pathways imported from KEGG.
In summary, Thanos provides functions to import metagenomic data from standard formats, perform a fine-grained functional annotation of the genes within the samples, and visualize the resulting profiles, potentially aggregated into an annotated reaction graph.
The functionality and usage of the package are described in the next sections.

# Methods (software description)

The overall design of the package draws inspiration from `phyloseq`.
Briefly, the main component of a `phyloseq` object is the OTU abundance table, in the form of a numeric matrix with taxa on the rows and samples on the columns (or vice versa).
The OTU table can be optionally decorated with sample metadata, an expanded taxonomy table, a phylogenetic tree, and even the reference genome of each taxon.
Naturally, `phyloseq` also provides functions to perform common manipulations on the abundance table, such as filtering samples, pruning taxa, or aggregating by taxonomy.
The main idea behind Thanos is that the same objects and methods can be used on *gene* abundances as well as taxa abundances.
Thus, the main purpose of our package is to perform a mapping from the taxonomy space to the functional space of a metagenomic sample, while using `phyloseq` objects to keep track of the abundances.
We therefore inherit all of the useful methods that have already been implemented in `phyloseq`.

We consider a sample to be a set of DNA sequences.
As Thanos doesn't require (but can still make use of) the taxonomy information, the sequences can be either assembled contigs or binned MAGs.
Each contig or MAG is associated with a number representing its sequencing depth in each sample, and is initially stored in a `phyloseq` object.
In order to map the sample to its functional space, Thanos uses the HMMER software to perform a sequence similarity search.
However, unlike existing tools like Prokka or the EggNOG-mapper, which search each gene in the sample against a database of target sequences, Thanos reverses this process, uses the genes of interest as query and the sample as target database.
The advantage of this approach is that it lets researchers build a custom profile hidden markov model (HMM) for their gene of interest, allowing them to capture the specific sequence variation that they are interested in.
Nevertheless, as custom profiles are not always necessary, Thanos also has the ability to automatically build a profile for a given gene by leveraging the KEGG orthologs database.

Once the profiles for the genes of interest have been generated, Thanos searches for them across the samples, keeping track of the contigs/MAGs where each gene is found, and summing their depths.
Next, Thanos performs another `HMMER` search, this time not for a gene of interest, but for a control gene: a universal marker gene that is conserved across all microorganisms (bacteria, archaea, or both, depending on the research question).
Dividing the total depth of a gene of interest by the total depth of the control gene within each sample produces a normalized score that can be interpreted as the average copy-number of the gene of interest in the sample, and can therefore be compared across independent samples.
If taxonomic information is available, the depth calculation can be done for each taxon independently.
The results are compiled into another `phyloseq` matrix, where instead of OTUs we have functional profiles.

Finally, Thanos provides functions to visualize the results using the popular `ggplot2` package.
There are three types of plot: bar plots, which show the depth profile across samples, potentially grouped by taxon; box plots, displaying aggregate statistics about groups of samples; and annotated reaction graphs, reproducing a KEGG module and coloring the reactions by the depths of the enzymes that catalyze them.
In the subsequent subsections, the implementation will be discussed in detail.

## Importing depths

First, we need to import the pre-computed depths.
These files can be obtained from standard tools in the metagenomic arsenal, such as MetaBAT or CoverM, and are generated as part of a standard workflow.
Thanos can read depths files in these standard formats.
In the case of MAG depths, this is typically a single file with all of the MAGs in the rows and the samples on the columns, and can be imported directly as a `phyloseq` object.
For contig depths, due to their large size, there is usually one file for each sample, where contigs from that sample are measured across all other samples; Thanos simply takes the list of files and concatenates them into a single `phyloseq` object.
One complication is that the header of the contig depths files often has a prefix or suffix denoting the sample where the contig comes from; since each file would have a different header, this prevents the concatenation.
For this reason the user has to specify a pattern and a replacement that will be applied to the headers, so that the sample-specific part can be removed.

## Deploying HMMER

Next, we need a profile HMM file for each gene of interest.
If the user has already produced or downloaded custom HMMs, nothing else is needed.
Otherwise, HMMs can be generated directly from KEGG orthologs: the user just needs to supply a KO identifier, and all the genes in the orthologous family are downloaded, aligned, and converted into an HMM profile.
The download uses the `KEGGREST` package, which interfaces with the KEGG API.
For the alignment we rely on the `msa` package from Bioconductor, which offers several choices for the algorithm, with "Muscle" being the default.
Once the alignment is obtained, it is converted into an HMM profile using `HMMER`.
No matter how it was generated, the HMM files for the genes of interest must then be searched into the samples.

For this purpose, the user is expected to provide sequence files in FASTA format containing the called genes from each entry in the previously generated `phyloseq` depths object.
In particular, when dealing with MAG depths, there should be one fasta file for each MAG.
Since each individual MAG normally has many genes, there must be a mechanism to associate each gene to the MAG where it comes from.
Thanos uses a linker function


## Aggregating and normalising

## Visualization

# Results

## Mags workflow

## Contigs workflow

# Discussion

An obvious caveat is that, even if a gene has a high DNA copy number, this does not necessarily mean that the gene will be highly expressed.
Thus, whenever possible, metagenomics data should always be complemented by meta-transcriptomic experiments.

# Overview of the method

thanos requires three inputs: a list of genes of interest, which are identified by their KO (KEGG Ortholog) number [@kanehisa2000kegg]; the depths file, which represents the abundances of OTUs across samples (either raw contigs depths or binned MAGs depths); and the sequences files, which associate each OTU with the protein sequences that it contains. The `thanos` method consists in looking for the gene of interest in the sequences file using HMMer [@eddy2011hmm], then mapping the results back to the OTUs, in order to get a depth profile of the gene of interest across samples. However, these raw depths are not comparable across samples due to different overall sequencing depths. For this reason, they are normalised by the depths of the hits of single-copy marker genes that are universally conserved (for instance, any of the 120 marker genes from GTDB [@parks2021gtdb]). This enables us to interpret the final score as the average copy number of the gene of interest in the sample. An overview of the functions available in the package is shown in \autoref{fig:workflow}. However, there are also high-level function that automate most of the analysis in one step, as well as functions optimised for large-scale data. The reader is invited to consult the vignette of the package for further details.

![thanos workflow.\label{fig:workflow}](workflow.png)

# Acknowledgments:
This work was supported by the China Scholarship Council and the Science & Technology Basic Resources Investigation Program of China (Grant No. 2017FY100300).
