---
date: 5 May 2024
bibliography: paper.bib
---

# Introduction

The field of metagenomics has experienced significant growth over recent decades [@zhang_2021], offering unprecedented insights into microbial communities across various environments, ranging from the human gut [@fzaal_2022] to marine ecosystems [@sunagawa_2020].
The progress was largely enabled by the advent of next-generation sequencing technologies and the development of algorithms for the reconstruction of individual microbial genomes from the pooled and fragmented DNA extracted from environmental samples [@zhang_2021; @yang_2021].
As is often the case, however, collecting more data does not necessarily lead to an improved understaing of the underlying biological systems.
Thus, with increasing sequencing data, there also increases the need for tools capable of interpreting them.

One of the most basic questions that can be asked concerns which taxa are present in a sample.
The `phyloseq` [@mcmurdie_2013] package is a powerful tool for the exploration of microbiome profiles, offering a practical approach to investigate the taxonomic composition of metagenomic samples collected from diverse environments.
However, the taxonomic profile offers only one view into the complex and multifaceded nature of biological samples.
The gene composition of a sample offers a complementary view, one that can help answer questions such as: does carbon fixation occur in this environment?, or: is methane metabolism more active in lakes or in the Atlantic ocean?
This sort of information is not always reflected in the taxonomic profiles [@franzosa_2018; yue_2023].
Furthermore, taxonomic profiling requires either reconstructing the full metagenome-assembled genomes (MAGs), which leads to a "binning bias" because unbinned sequences are not analyzed, or at least recovering several marker genes from the same operational taxonomic unit (OTU) [@milanese_2019], which can be time-consuming.

As the functional profiling task is somewhat more complex than taxonomic profiling, there is currently no standard way to perform this analysis.
Tools like Prokka [@seeman_2014] or the EggNog-mapper [@cantalapiedra_2021] can indeed provide a bulk-level overview of the functional composition by performing sequence-similarity searches for each gene in the sample, but researchers are often interested in more specific questions about individual genes or metabolic pathways, and require a more in-depth analysis [@yue_2023].

In this work, we present Thanos, an R package that offers a convenient way to perform functional profiling with a gene- or pathway- centric approach.
The package provides quantitative functional information through a "depth score" for each gene of interest across samples: genes with a higher depth are thought to be more prevalent and active in the sample.
Moreover, we introduce a normalization strategy that makes depth scores comparable across samples and even across independent sequencing projects, enabling the comparative analysis of multiple environments at the same time.
The depth scores of individual genes can also be aggregated into their natural higher-level units, the metabolic pathways imported from KEGG.

# Methods

As R is one of the most popular languages for bioinformatics, not least thanks to the Bioconductor project, Thanos was implemented in R and makes use of its ecosystem.
The overall design of the package draws inspiration from `phyloseq` [@mcmurdie_2013].
Briefly, the main component of a `phyloseq` object is the OTU abundance table, in the form of a numeric matrix with taxa on the rows and samples on the columns (or vice versa).
The OTU table can be optionally decorated with sample metadata, an expanded taxonomy table, a phylogenetic tree, and even the reference genome of each taxon.
Naturally, `phyloseq` also provides functions to perform common manipulations on the abundance table, such as filtering samples, pruning taxa, or aggregating by taxonomy.
The main idea behind Thanos is that the same objects and methods can be used on *gene* abundances as well as taxa abundances.
Thus, the main purpose of our package is to perform a mapping from the taxonomy space to the functional space of a metagenomic sample, while using `phyloseq` objects to keep track of the abundances.
We therefore inherit all of the useful methods that have already been implemented in `phyloseq`.

We consider a sample to be a set of DNA sequences.
As Thanos doesn't require (but can still make use of) the taxonomy information, the sequences can be either assembled contigs or binned MAGs.
Each contig or MAG is associated with a number representing its sequencing depth in each sample, and is initially stored in a `phyloseq` object.
In order to map the sample to its functional space, Thanos uses the `HMMER` software [@eddy_2011] to perform a sequence similarity search.
However, unlike existing tools like Prokka [@seeman_2014] or the EggNOG-mapper [@cantalapiedra_2021], which search each gene in the sample against a database of target sequences, Thanos reverses this process, using the genes of interest as query and the sample as target database.
The advantage of this approach is that it lets researchers build a custom profile hidden markov model (HMM) for their gene of interest, allowing them to capture the specific sequence variation that they are interested in.
Nevertheless, as custom profiles are not always necessary, Thanos also has the ability to automatically build a profile for a given gene by leveraging the KEGG orthologs database [@kaneisha_2000].

Once the profiles for the genes of interest have been generated, Thanos searches for them across the samples, keeping track of the contigs/MAGs where each gene is found, and summing their depths.
Next, Thanos performs another `HMMER` search, this time not for a gene of interest, but for a control gene: a universal marker gene that is conserved across all microorganisms (bacteria, archaea, or both, depending on the research question).
Dividing the total depth of a gene of interest by the total depth of the control gene within each sample produces a normalized score that can be interpreted as the average copy-number of the gene of interest in the sample, and can therefore be compared across independent samples.
If taxonomic information is available, the depth calculation can be done for each taxon independently.
The results are compiled into another `phyloseq` matrix, where instead of OTUs we have the gene or pathway of interest.

Finally, Thanos provides functions to visualize the results using the popular `ggplot2` package [@wickham_2016].
There are three types of plot: bar plots, which show the depth profile across samples, potentially grouped by taxon; boxplots, displaying aggregate statistics about groups of samples; and annotated reaction graphs, reproducing a KEGG module and coloring the reactions by the depths of the enzymes that catalyze them [@kaneisha_2000].
Figure @fig:workflow gives a global overview of the packages' functionality, while in the subsequent subsections the implementation will be discussed in detail.

![Overview of the Thanos workflow.](../workflow.png){#fig:workflow width=100%}

## Importing depths files

First, we need to import the pre-computed depths.
These files can be obtained from standard tools in the metagenomic arsenal, such as MetaBAT [@kang_2019] or CoverM [@coverm_github], and are generated as part of a standard workflow [@yang_2021].
Thanos can read depths files in these standard formats.
In the case of MAG depths, this is typically a single file with all of the MAGs in the rows and the samples on the columns, and can be imported directly as a `phyloseq` object.
For contig depths, due to their large size, there is usually one file for each sample, where contigs from that sample are measured across all other samples; Thanos simply takes the list of files and concatenates them into a single `phyloseq` object.
One complication is that the header of the contig depths files often has a prefix or suffix denoting the sample where the contig comes from; since each file would have a different header, this prevents the concatenation.
For this reason the user has to specify a pattern and a replacement that will be applied to the headers, so that the sample-specific part can be removed.

## Running HMMER

Next, we need a profile HMM file for each gene of interest.
Typcal genes of interest are enzymes that perform key metabolic functions.
If the user has already produced or downloaded custom HMMs, nothing else is needed.
Otherwise, HMMs can be generated directly from KEGG orthologs: the user just needs to supply a KO identifier, and all the genes in the orthologous family are downloaded, aligned, and converted into an HMM profile.
The download uses the `KEGGREST` Bioconductor package, which interfaces with the KEGG API.
For the alignment we rely on the `msa` package from Bioconductor, which offers several choices for the algorithm, with "Muscle" being the default.
Once the alignment is obtained, it is converted into an HMM profile using the local `HMMER` binaries.

No matter how it was generated, the HMM files for the genes of interest must then be searched into the samples.
For this purpose, the user is expected to provide sequence files in FASTA format containing the called genes from each entry in the previously generated `phyloseq` depths object.
In particular, when dealing with MAG depths there should be one fasta file for each MAG, and when dealing with contig depths there should be one fasta file for each sample.
These files are produced by popular tools like `Prodigal` [@hyatt_201] or `Prokka` [@seeman_2014], which are already included in comprehensive metagenomic pipelines like `nf-core/mag` [@krakau_2022; ewels_2020].
Thanos provides the function `search_hmm(hmm_file, target_fasta_files)`, with which spawns an `hmmsearch` process and parses the results.
In this case, the results will be the IDs of the genes that bear sequence similarity with the HMM profile, along with their respective scores.
Users can specify the minimum score threshold to retain hits.
The search should be repeated for the control gene, for which Thanos already provides an HMM file, but users can also use their own if they wish.
The default control gene is GrpE, a nucleotide exchange factor that is important for protein folding and heat-shock response [@bracher_2015; delaney_1990].
The HMM profile for this gene was derived from GTDB v214 [@chaumeil_2022; @parks_2022] marker files: `bac120_r214_reps_PF01025.20.afa`.

## Aggregation and normalisation

At the end of the search phase we thus have a list of genes that are homologous to the given HMM profile, as well as a list of genes that are homologous to the control HMM profile.
Since each sample usually contains many MAGs or contigs, and each individual MAG or contig contains many genes, there must be a mechanism to associate each gene to the MAG or contig where it comes from.
Thanos uses a linker function, which can also be specified by the user, to achieve this.
The linker function takes the name of a FASTA file and the ID of the gene, and returns the MAG or contig where it comes from.
Through the linker function it becomes possible to filter the MAGs or contigs according to whether they contain the target gene.
For each sample, Thanos aggregates the depths of all the MAGs or contigs that contain the gene of interest, and divides it by the aggregated depth of all the MAGs or contigs that contain the control gene.
The resulting score represent how prevalent is the gene of interest compared to a universal single copy gene, in each sample.
For the convenience of the user, two linker functions which cover the most common cases are already built-in.
As an additional feature, when the taxonomy assignments of the MAGs are known, it is possible to stratify the score by taxonomy.
This means that there will be one score for each taxon in each sample, making it easy to make hypotheses about the role of a particular taxon in the ecosystem.
All the scores are saved in a `phyloseq` object, which is returned to the user.

## Visualization

The last part of the Thanos workflow consists in the visualization of the results.
Advanced users can of course compose their own plots starting from the aggregated results, but three plot types are also provided by default.
The first is a bar plot of the depth scores.
We provide a flexible interface where users can choose what to show on the *x*-axis, and the depths will be automatically aggregated by that variable.
Indeed, as the results are normal `phyloseq` objects, they can be decorated with sample metadata or taxonomy tables.
By default, samples are on the *x*-axis, but users may choose to aggregate the samples into subgroups, or to show the depths by taxonomy instead.
The second plot type is a box plot, useful to compare groups of samples.
Again, users have all the freedom to customize the groups.
Finally, we provide a function that plots the KEGG reaction graph of a whole module, where each enzyme is colored by its depth score in the samples of interest.

## Parallelization

As metagenomic datasets can be rather big, Thanos makes it possible to run the HMM searches in parallel, which can dramatically speed-up the code execution.
There are two nested levels of parallelism: firstly, users can control how many parallel `hmmsearch` processes are spawned, and secondly, for each process it is possible to choose how many threads it will use.
The outermost level of parallelism is most useful when there are many protein sequence databases, whereas the innermost level is especially useful when the individual sequence databases are large.
Reading contigs depths files is also parallelized.

## Dependencies

Thanos depends on R and the following packages: `phyloseq`, `data.table`, `KEGGREST`, `msa`, `Biostrings`, and `ggplot2`.
In addition, the `HMMER` binaries must be installed separately.

# Results

To illustrate the functionality of our software, we will showcase two applications.
The data come from the TARA ocean project [@sunagawa_2015] (Suppl. Table 1) and they were pre-processed as follows.
First, we downloaded the raw reads for two ocean "provinces" (according to the nomenclature of the original publication): the Red Sea and the Mediterranean Sea, corresponding to 11 distinct sampling stations.
Then we used the `nf-core/mag` v2.5.4 automated pipeline for assembly, binning, and annotation [@krakau_2022].
This provided all the files necessary for running Thanos.

## Mags workflow: sulfur metabolism by taxonomy

Our aim was to investigate which taxa have the potential to perform sulfur metabolism in these two seas.
Thanos minimally requires three inputs: a list of genes of interest, the MAG depth files, and the protein sequence files, all generated by the `nf-core/mag` workflow.
Since we wanted to stratify the analysis by taxonomy, we also provided a table with the GTDB taxonomy of each MAG, also generated by `nf-core/mag`.
As for the genes, we extracted all the genes in KEGG's "Assimilatory sulfate reduction" pathway.
We just had to give Thanos the KEGG ortholog IDs of the genes and the paths to the files generated during the `nf-core/mag` workflow.

First, we can explore the sulfur assimilation potential by plotting the depth score of the genes across stations (Figure @fig:mags A).
This give us a feeling for any overall differences between the stations, as well as which are the most important phyla.
In this case, the abundances of sulfur genes is relatively uniform across samples, except possibly for station `TARA_022`, which shows lower abundances.
Moreover, the contribution of phyla to sulfur assimilation also appear uniform, with no prevalent taxon.
The most abundant enzymes are APR, CysC, CysN, and CysCN.

If we are especially interested in one enzyme, say, CysCN, we can explore it more in-depth by comparing its abundance across two environments.
We found that its abundance is significantly lower in the Red sea than in the Mediterranean sea (Figure @fig:mags B).
Nevertheless, the relative contributions of the taxa are not greatly different, with the exception of station `TARA_022` (Figure @fig:mags C).

![Bla](../mags_patchwork.png){#fig:mags width=100%}

## Contigs workflow: prevalence of glycolysis

In the second example, we examine the prevalence of glycolisys genes.
As we are not interested in the taxonomy, we can use the contigs rather than the MAGs, so that we retain even the unbinned DNA.
We can annotate the reaction graph of glycolysis with the depths scores computed for our samples (Figure @fig:contigs A).
We see that two reactions, namely glucose to glucose-6-phosphate performed by the glucose phosphotransferase enzyme, and glyceraldeide-3-phosphate to 3-phosphoglycerate performed by glyceraldehyde-3-phosphate dehydrogenase with ferredoxin cofactor, are almost absent.
On the other hand, the enzymes phosphoglycerate kinase, catalyzing the reaction from 3-phosphoglyceroyl-phosphate to to 3-phosphoglycerate, have an average copy number of 2.8 in these samples.
As this reaction can be performed by four different enzymes, we also investigated the abundance of each individual gene (Figure @fig:contigs B), finding that GapB (K00150) and GapA (K00134) are the most abundant, whereas gapor (K11389) is virtually inexistent.

![Bla](../contigs_patchwork.png){#fig:contigs width=100%}

# Discussion

We developed a package to streamline a pathway-centric analysis of metagenomics data.
It can analyze both contig-level data and MAG-level data within a single, general framework.
The software is user-friendly and efficient, and it integrates well within the existing R ecosystem, in particular with the `phyloseq` and `ggplot2` packages.
An obvious caveat of the Thanos method is that, even if a gene has a high DNA copy number, this does not necessarily mean that the gene will be highly expressed.
Thus, whenever possible, metagenomics data should be complemented by meta-transcriptomic or even proteomics experiments.
Even with this limitation, functional profiling of metagenomic samples can offer insights into not just which bacteria populate an environment, but also what are they doing there and how much are they doing it.
This is especially relevant as even different strains of the same species can harbor sometimes vastly different gene portfolios, and therfore perform vastly different metabolic reactions [@kaper_2004; brockhurst_2019].
Thus, by looking at the genes directly instead of the taxonomy, we can gain a more nuanced picture of the role of the microbial community in an environment.
Indeed, our package has been optimized for comparing functional profiles across environments.
Even though we focused on metabolic pathways in this paper, the applicability of Thanos extends to any other phenotype that can be linked to genes.
For example, it is possible to compare the abundance of pathogenic genes, or antimicrobial genes, or even different classes of CRISPR systems.
In conclusion, we hope that our framework will enable new discoversies in environmental research.

# Code availability

The Thanos package and the code used to produce the plots in this paper are hosted on GitHub at the url [https://github.com/zhezhaozoe/thanos].

# Acknowledgments:
This work was supported by the China Scholarship Council and the Science & Technology Basic Resources Investigation Program of China (Grant No. 2017FY100300).
