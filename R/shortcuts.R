#' Perform HMM Search on a KEGG Orthology (KO) Alignment
#'
#' Given a KO identifier, this function retrieves the corresponding protein sequences and builds a Multiple Sequence Alignment (MSA), then builds a Hidden Markov Model (HMM) based on the MSA, and then searches
#' the HMM against a set of sequence databases. The function allows adjusting the search method,
#' computational resources, and inclusion threshold.
#'
#' @param ko A character string representing the KEGG Orthology (KO) identifier.
#' @param method Character. The method used for performing the MSA. Defaults to "Muscle".
#' @param cpu Integer. The number of CPUs to use for the search. Defaults to 1.
#' @param incE Numeric. The inclusion threshold for the HMM search. Defaults to 1e-6.
#' @param ... Additional arguments passed to kegg_kegg_msa().
#'
#' @return The result of the HMM search.
#' @examples
#' tblout_from_ko("K00001")
#' @export
tblout_from_ko <- function(ko, dbs, method = "Muscle", cpu = 1, incE = 1e-6, ...) {
  aln <- get_kegg_msa(ko, method = method, ...)
  hmm <- build_hmm(aln)
  search_hmm(hmm, dbs, cpu = cpu, incE = incE)
}

#' Perform HMM-based search from an alignment
#'
#' This function takes aligned fasta sequences (AFA) as an input, constructs
#' a Hidden Markov Model (HMM), and uses it to search in designated databases.
#'
#' @param afa Aligned fasta sequences in string format or an object that
#'  can be interpreted as such.
#' @param cpu The number of CPU cores to use for the search process.
#'  Defaults to 1 to ensure compatibility with all systems.
#' @param incE The inclusion threshold E-value for reported hits in the
#'  search. Lower values are more stringent. Defaults to 1e-6.
#' @param ... Additional arguments to be passed to the search function.
#'
#' @return Returns the path to the tblout results.
#'
#' @examples
#' # Example assuming 'sequences.afa' is your aligned fasta file
#' afa_content <- readLines("sequences.afa")
#' search_results <- tblout_from_afa(afa_content)
#' print(search_results)
#'
#' @export
tblout_from_afa <- function(afa, dbs, cpu = 1, incE = 1e-6) {
  hmm <- build_hmm(afa)
  search_hmm(hmm, dbs, cpu = cpu, incE = incE)
}

#' Perform a Multiple Sequence Alignment and Find Hits in Databases
#'
#' This function first reads a fasta file containing amino acid sequences, performs
#' a multiple sequence alignment using the specified method, and then searches
#' for hits in databases using the aligned sequences. The results are returned
#' in a specific format.
#'
#' @param faa Character. The path to the fasta file containing amino acid sequences.
#' @param dbs Character or Character vector. Specifies the databases against which the
#' sequences will be searched.
#' @param method Character. The default method for sequence alignment is "Muscle".
#' Other methods supported by `msa` function can be used.
#' @param cpu Integer. The number of CPU threads to use for the search. The default is 1.
#' @param incE Numeric. The inclusion threshold E-value for considering a database hit
#' significant. Defaults to 1e-6.
#' @param ... Additional arguments passed to the `msa` function.
#'
#' @return Returns the result of searching the aligned sequences against specified
#' databases. The format of the returned object depends on the implementation of
#' `tblout_from_afa`.
#'
#' @examples
#' \dontrun{
#' tblout_from_faa("sequences.faa", dbs=c("db1", "db2"), method="ClustalW")
#' }
#'
#' @export
#' @importFrom Biostrings readAAStringSet
#' @importFrom msa msa
tblout_from_faa <- function(faa, dbs, method = "Muscle", cpu = 1, incE = 1e-6, ...) {
  mySequences <- Biostrings::readAAStringSet(faa)
  aln <- msa::msa(mySequences, method = method, ...)
  tblout_from_afa(aln, cpu = cpu, incE = incE)
}

#' Obtain Contigs Hits Depths from Depth Files
#'
#' This function processes depth files to extract and calculate the hits depths for contigs.
#' Contig depths files are typically large (several hundred MB per sample), therefore, with many samples, it's impractical to load all of the data at once in memory. This function provides a way to calculate the depths sample-wise, to avoid storing all the data in memory.
#'
#' @param depths_files A character vector of file paths to the depths files.
#' @param pattern Character. A regular expression pattern to match in contigs names within the depth files.
#' @param replacement Character. A string to replace the matched pattern in contig names.
#' @param query_tblout File path to the tblout file containing query hits.
#' @param control_tblout File path to the tblout file containing control hits.
#' @param linker A function to link contig identifiers between the query and control tblout files.
#' @param verbose Logical. Indicates whether to print detailed processing messages. Defaults to FALSE.
#' @return A phyloseq object with the relative depths of query vs control genes across all `detphs_files`.
#' @examples
#' get_contigs_hits_depths(c("depths1.txt", "depths2.txt"), "_contig", "_sequence", "query.tblout", "control.tblout")
#' @export
get_contigs_hits_depths <- function(depths_files, pattern, replacement, query_tblout, control_tblout, linker = contigs_linker, verbose = F) {
  # We use seq_along here because iterating over the depth_files directly
  # the name attribute is lost, and the name is used in import_contig_depths
  mers <- lapply(seq_along(depths_files), function(i) {
    cdf <- depths_files[i]
    if (isTRUE(verbose)) {
      message(cdf)
    }
    otus <- import_contig_depths(cdf, pattern, replacement)
    get_hits_depths(otus, query_tblout, control_tblout, linker, phyloseq = F)
  })
  mer <- rbindlist(mers, use.names = T)[, lapply(.SD, sum), by = "ID"]
  res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
  setnames(res, names(res), sub("\\.x$", "", names(res)))
  m <- as.matrix(res, rownames.value = mer$ID)
  # XXX: set NaNs (0/0) to zero
  m[is.nan(m)] <- 0
  otu_table(m, taxa_are_rows = T)
}

#' Retrieve depths of hits from KEGG Orthology entries across samples
#'
#' This function takes KEGG Orthology (KO) identifiers, a phyloseq object, a set of databases,
#' a control tblout formatted file, and computes
#' the depths of hits related to the input KOs in the provided sample set.
#'
#' @param kos A character vector of KEGG Orthology (KO) identifiers.
#' @param ps A phyloseq object containing the samples data.
#' @param dbs A character vector indicating the databases to be used.
#' @param control_tblout A tblout formatted file with the control genes' hits.
#' @param linker Function that links OTU identifiers between tblout tables and the phyloseq object. The function must take two arguments: SeqFile and Target.
#' @param msa_method The method used for multiple sequence alignment. "Muscle" by default.
#' @param cpu The number of CPUs to be used for computation. Integer value. 1 by default.
#' @param incE The inclusion E-value threshold for considering a hit significant. Numeric value. 1e-6 by default.
#' @param taxrank Optional; character specifying the taxonomic rank at which to aggregate hits. Valid options depend on the taxonomic ranks present in the `ps` object. If not specified, no aggregation is performed.
#'
#' @return A list where each element corresponds to a KO identifier inputted and contains
#'         the computed depths of hits for that KO in the sample set.
#'
#' @examples
#' kos <- c("K00001", "K00002")
#' ps <- phyloseq_object # a pre-loaded or created phyloseq object
#' dbs <- c("db1", "db2")
#' control_tblout <- "control_hits.tblout"
#' linker <- read.csv("linker_table.csv")
#' hit_depths <- get_hits_depths_from_kos(kos, ps, dbs, control_tblout, linker)
#'
#' @export
get_hits_depths_from_kos <- function(kos, ps, dbs, control_tblout, linker, msa_method = "Muscle", cpu = 1, incE = 1e-6, taxrank = NULL, ...) {
  if (is.null(names(kos))) {
    names(kos) <- kos
  }
  sapply(simplify = FALSE, kos, function(ko) {
    query_tblout <- tblout_from_ko(ko, dbs, method = msa_method, cpu = cpu, incE = incE, ...)
    get_hits_depths(ps, query_tblout, control_tblout, linker, taxrank = taxrank)
  })
}

#' Get Contigs Hits Depths from KEGG Orthology (KO) Identifiers
#'
#' For each KEGG Orthology (KO) identifier provided, this function retrieves contig hits and their depths from the
#' specified depths files.
#'
#' @param kos A character vector of KEGG Orthology (KO) identifiers.
#' @param depths_files A character vector indicating the paths to depths files.
#' @param pattern A regex pattern used to match and extract contig names from the depths files.
#' @param replacement A character string to replace in the extracted contig names, based on the `pattern` argument.
#' @param dbs A not used parameter in the current implementation.
#' @param control_tblout Character string indicating the path to the tblout file from the control.
#' @param linker A function that links contig names to their hits and depths. Defaults to `contigs_linker`.
#' @param msa_method Character string representing the method used for multiple sequence alignment. Defaults to `"Muscle"`.
#' @param cpu An integer indicating the number of CPUs to be used in computations. Defaults to `1`.
#' @param incE A numeric value specifying the inclusion E-value threshold. Defaults to `1e-6`.
#' @param verbose A logical indicating whether to display detailed output. Defaults to `FALSE`.
#'
#' @return Returns a list where each element corresponds to the result from a single KO identifier, containing the contig hits and their depths information.
#'
#' @examples
#' \dontrun{
#' kos <- c("K00001", "K00002")
#' depths_files <- c("sample1.depths", "sample2.depths")
#' pattern <- "_(\d+)$"
#' replacement <- ""
#' control_tblout <- "control.tblout"
#' results <- get_contigs_hits_depths_from_kos(kos, depths_files, pattern, replacement, NULL, control_tblout)
#' }
#'
#' @export
get_contigs_hits_depths_from_kos <- function(kos, depths_files, pattern, replacement, dbs, control_tblout, linker = contigs_linker, msa_method = "Muscle", cpu = 1, incE = 1e-6, verbose = F) {
  if (is.null(names(kos))) {
    names(kos) <- kos
  }
  sapply(simplify = FALSE, kos, function(ko) {
    query_tblout <- tblout_from_ko(ko, method = msa_method, cpu = cpu, incE = incE)
    get_contigs_hits_depths(depths_files, pattern, replacement, query_tblout, control_tblout, linker = linker, verbose = verbose)
  })
}
