#' Build an HMM profile for the given KO identifiers
#'
#' Given a KO identifier, this function retrieves the corresponding protein sequences and builds a Multiple Sequence Alignment (MSA), then builds a Hidden Markov Model (HMM) based on the MSA.
#'
#' @param ko A character string representing the KEGG Orthology (KO) identifier.
#' @param method Character. The method used for performing the MSA. Defaults to "Muscle".
#' @param ... Additional arguments passed to kegg_kegg_msa().
#'
#' @return The path to the HMM file.
#' @examples
#' hmm_from_ko("K00001")
#' @export
build_hmm_from_ko <- function(kos, method = "Muscle", hmmer_path = "", ...) {
  if (is.null(names(kos))) {
    names(ko) <- kos
  }
  sapply(simplify = FALSE, kos, function(ko) {
    build_hmm(get_kegg_msa(ko, method = method, ...), hmmer_path = hmmer_path)
  })
}

#' Obtain Contigs Hits Depths from Depth Files
#'
#' This function processes depth files to extract and calculate the hits depths for contigs.
#' Contig depths files are typically large (several hundred MB per sample), therefore, with many samples, it's impractical to load all of the data at once in memory. This function provides a way to calculate the depths sample-wise, to avoid storing all the data in memory.
#'
#' @param depths_files A character vector of file paths to the depths files.
#' @param sub_pattern Character. A regular expression pattern to match in contigs names within the depth files.
#' @param sub_replacement Character. A string to replace the matched pattern in contig names.
#' @param query_tblout File path to the tblout file containing query hits.
#' @param control_tblout File path to the tblout file containing control hits.
#' @param linker A function to link contig identifiers between the query and control tblout files.
#' @param verbose Logical. Indicates whether to print detailed processing messages. Defaults to FALSE.
#' @return A phyloseq object with the relative depths of query vs control genes across all `detphs_files`.
#' @examples
#' get_contigs_hits_depths(c("depths1.txt", "depths2.txt"), "_contig", "_sequence", "query.tblout", "control.tblout")
#' @export
get_contigs_hits_depths <- function(depths_files,
sub_pattern, sub_replacement, query_tblout, control_tblout,
linker = contigs_linker, verbose = F) {
  # We use seq_along here because iterating over the depth_files directly
  # the name attribute is lost, and the name is used in import_contig_depths
  mers <- lapply(seq_along(depths_files), function(i) {
    cdf <- depths_files[i]
    if (isTRUE(verbose)) {
      message(cdf)
    }
    otus <- import_contig_depths(cdf, sub_pattern, sub_replacement)
    get_hits_depths(
      otus,
      query_tblout, control_tblout,
      linker,
      phyloseq = F)
  })
  mer <- rbindlist(mers, use.names = TRUE)[, lapply(.SD, sum), by = "ID"]
  merged_hits_to_otu_table(mer, taxa_are_rows = TRUE)
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
get_hits_depths_from_kos <- function(kos, ps, dbs,
control_tblout, linker, msa_method = "Muscle",
cpu = 1, incE = 1e-6, taxrank = NULL, ...) {
  .Deprecated("get_hits_depths_from_hmm()", package = "thanos")
  if (is.null(names(kos))) {
    names(kos) <- kos
  }
  sapply(simplify = FALSE, kos, function(ko) {
    query_tblout <- tblout_from_ko(
      ko,
      dbs,
      method = msa_method,
      cpu = cpu,
      incE = incE,
      ...)
    get_hits_depths(
      ps,
      query_tblout, control_tblout,
      linker,
      taxrank = taxrank)
  })
}

#' Retrieve depths of hits from KEGG Orthology entries across samples
#'
#' @description This function executes a HMMER search against a database for a set of query and control HMM profiles. It returns the hits depths from KEGG Orthology (KO) identifiers for each query relative to the control.
#'
#' @param queries_hmm A character vector of paths to the query HMM files.
#' @param control_hmm A character string specifying the path to the control HMM file.
#' @param ps A parameter set or data required by `get_hits_depths` function.
#' @param dbs A character vector specifying the paths to the database files to be searched.
#' @param linker A parameter or function used to connect or process the information between the query and control results.
#' @param cpu The number of CPUs to use for the HMMER search. Defaults to 1.
#' @param incE The inclusion E-value threshold for the HMMER search. Defaults to 1e-6 to filter significant hits.
#' @param taxrank An optional parameter specifying the taxonomic rank consideration. Can be used within `get_hits_depths` function.
#' @param ... Additional arguments passed to the `get_hits_depths` function.
#'
#' @return A list where each element corresponds to the result of `get_hits_depths` for a query HMM, containing information on the hits' depths based on KO identifiers.
#'
#' @examples
#' # Assuming proper environment setup and availability of required files
#' queries = c("path/to/query1.hmm", "path/to/query2.hmm")
#' control = "path/to/control.hmm"
#' db_paths = c("path/to/db1", "path/to/db2")
#' # Example call
#' results = get_hits_depths_from_hmm(queries, control, parameter_set, db_paths, linker_function)
#'
#' @export
get_hits_depths_from_hmm <- function(queries_hmm, control_hmm, ps, dbs, linker, taxrank = NULL, parallel_processes = 1, cpus_per_process = 1, incE = 1e-6, hmmer_path = "") {
  control_tblout <- search_hmm(control_hmm, dbs, parallel_processes = parallel_processes, cpus_per_process = cpus_per_process, incE = incE, hmmer_path = hmmer_path)
  sapply(simplify = FALSE, queries_hmm, function(query_hmm) {
    query_tblout <- search_hmm(query_hmm, dbs, parallel_processes = parallel_processes, cpus_per_process = cpus_per_process, incE = incE, hmmer_path = hmmer_path)
    get_hits_depths(
      ps,
      query_tblout, control_tblout,
      linker,
      taxrank = taxrank)
  })
}

#' Get Contigs Hits and Depths from HMM Searches
#'
#' This function performs HMM searches for query and control sequences against a set of databases, retrieves hits, and calculates depths for contigs.
#' 
#' @param queries_hmm Character vector specifying the file paths to the query HMM profiles.
#' @param control_hmm Character string specifying the file path to the control HMM profile.
#' @param depths_files Character vector with paths to the depth files.
#' @param sub_pattern Character string representing the pattern to identify contig names within depth files.
#' @param sub_replacement Character string with the replacement for the sub_pattern in contig names.
#' @param dbs Character vector specifying the paths to the databases to search against.
#' @param linker Function that links hits to contigs. If not supplied, a default contigs_linker function will be used.
#' @param cpu Integer specifying the number of CPUs to use for HMM searches. Defaults to 1.
#' @param incE Numeric, inclusion E-value threshold for HMM searches. Defaults to 1e-6.
#' @param hmmer_path Character string specifying the path to the HMMER binaries. If empty, system's PATH will be used.
#' @param verbose Logical indicating whether to print detailed progress messages. Defaults to FALSE.
#'
#' @return A list where each element corresponds to the results from each query HMM, containing detailed hits and depths information for contigs.
#'
#' @examples
#' # Assuming appropriate files and HMM profiles exist:
#' result <- get_contigs_hits_depths_from_hmm(queries_hmm = c("/path/to/query1.hmm", "/path/to/query2.hmm"),
#'                                            control_hmm = "/path/to/control.hmm",
#'                                            depths_files = c("/path/to/depth1.tsv", "/path/to/depth2.tsv"),
#'                                            sub_pattern = "_[0-9]+$",
#'                                            sub_replacement = "",
#'                                            dbs = c("/path/to/db1", "/path/to/db2"),
#'                                            linker = custom_linker_function,
#'                                            cpu = 4, incE = 1e-10,
#'                                            hmmer_path = "/opt/hmmer/bin/",
#'                                            verbose = TRUE)
#'
#' @export
get_contigs_hits_depths_from_hmm <- function(
    queries_hmm, control_hmm,
    depths_files, sub_pattern, sub_replacement,
    dbs, linker = contigs_linker,
    parallel_processes = 1, cpus_per_process = 1, incE = 1e-6, hmmer_path = "",
    verbose = FALSE) {
  if (is.null(names(kos))) {
    names(kos) <- kos
  }
  control_tblout <- search_hmm(
    control_hmm,
    dbs,
    parallel_processes = parallel_processes,
    cpus_per_process = cpus_per_process,
    incE = incE,
    hmmer_path = hmmer_path
  )
  sapply(simplify = FALSE, kos, function(ko) {
    query_tblout <- search_hmm(
      queries_hmm,
      dbs,
      parallel_processes = parallel_processes,
      cpus_per_process = cpus_per_process,
      incE = incE,
      hmmer_path = hmmer_path
    )
    get_contigs_hits_depths(
      depths_files,
      sub_pattern,
      sub_replacement,
      query_tblout,
      control_tblout,
      linker = linker,
      verbose = verbose
    )
  })
}

#' Get Contigs Hits Depths from KEGG Orthology (KO) Identifiers
#'
#' For each KEGG Orthology (KO) identifier provided, this function retrieves contig hits and their depths from the
#' specified depths files.
#'
#' @param kos A character vector of KEGG Orthology (KO) identifiers.
#' @param depths_files A character vector indicating the paths to depths files.
#' @param sub_pattern A regex pattern used to match and extract contig names from the depths files.
#' @param sub_replacement A character string to replace in the extracted contig names, based on the `sub_pattern` argument.
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
#' sub_pattern <- "_(\d+)$"
#' sub_replacement <- ""
#' control_tblout <- "control.tblout"
#' results <- get_contigs_hits_depths_from_kos(kos, depths_files, sub_pattern, sub_replacement, NULL, control_tblout)
#' }
#'
#' @export
get_contigs_hits_depths_from_kos <- function(kos, depths_files,
sub_pattern, sub_replacement, dbs,
control_tblout, linker = contigs_linker,
msa_method = "Muscle", cpu = 1, incE = 1e-6, verbose = F, ...) {
  .Deprecated("get_contigs_hits_depths_from_hmm()", "thanos")
  if (is.null(names(kos))) {
    names(kos) <- kos
  }
  sapply(simplify = FALSE, kos, function(ko) {
    query_tblout <- tblout_from_ko(ko, dbs, method = msa_method, cpu = cpu, incE = incE, ...)
    get_contigs_hits_depths(depths_files, sub_pattern, sub_replacement, query_tblout, control_tblout, linker = linker, verbose = verbose)
  })
}
