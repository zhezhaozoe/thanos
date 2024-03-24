#' Build a Hidden Markov Model (HMM) from a sequence alignment
#'
#' This function takes a sequence alignment, converts it to the Stockholm format, and then uses the `hmmbuild` command
#' from the HMMER suite to build a Hidden Markov Model (HMM).
#'
#' @param aln A sequence alignment object. The alignment object should have an optional "name" attribute, which will be used as the name of the HMM.
#' @param hmmer_path Base path for HMMer binaries (if not in PATH).
#'
#' @return The path to the file containing the generated HMM.
#'
#' @examples
#' # Assuming `my_alignment` is an alignment object
#' hmm_file <- build_hmm(my_alignment)
#' # hmm_file now contains the path to the HMM file
#'
#' @note This function requires the external programs `esl-reformat` and `hmmbuild` from the HMMER suite to be
#'       installed and accessible in the system's PATH.
#'
#' @references
#' Eddy, S.R. "HMMER: biosequence analysis using profile hidden Markov models."
#' http://hmmer.org/
#'
#' @importFrom Biostrings writeXStringSet unmasked
#' @export
build_hmm <- function(aln, hmmer_path = "") {
  if (is.character(aln)) {
    if (!file.exists(aln)) {
      stop("File doesn't exist. Please provide a valid path")
    }
    afa <- aln
  } else if ("MsaAAMultipleAlignment" %in% class(aln)) {
    afa <- tempfile(fileext = ".afa")
    Biostrings::writeXStringSet(Biostrings::unmasked(aln), file = afa)
  } else {
    stop("Unrecognised format. aln can be either the path to an alignment in fasta format or an object of class from the msa package")
  }
  sto <- tempfile(fileext = ".sto")
  hmm <- tempfile(fileext = ".hmm")
  extargs <- if (!is.null(attr(aln, "name"))) {
    c("-n", attr(aln, "name"))
  } else {
    c()
  }
  esl_reformat_cmd <- if (nzchar(hmmer_path)) file.path(hmmer_path, "esl-reformat") else "esl-reformat"
  hmmbuild_cmd <- if (nzchar(hmmer_path)) file.path(hmmer_path, "hmmbuild") else "hmmbuild"
  system2(esl_reformat_cmd, c("stockholm", afa), stdout = sto)
  system2(hmmbuild_cmd, c(extargs, hmm, sto), stdout = NULL)
  return(hmm)
}

#' Perform HMMER search across multiple databases
#'
#' This function executes HMMER search on a given Hidden Markov Model (HMM) file across
#' multiple database files (dbs). It allows the setting of the number of CPUs to be used
#' and the inclusion E-value threshold. The results from searching each database are combined
#' into a single data table.
#'
#' @param hmm A string specifying the path to the HMM file to be searched with.
#' @param dbs A character vector, where each element is a path to a database file to be searched.
#' @param cpu An integer indicating the number of CPUs to be used for the search. Defaults to 1.
#' @param incE A numeric value specifying the inclusion E-value threshold. Defaults to 1e-6.
#' @param hmmer_path Base path for HMMer binaries (if not in PATH).
#'
#' @return A data frame containing the combined results of the HMM searches on all databases.
#'         Each row represents one HMMER hit, and the data frame includes a column 'SeqFile' indicating
#'         the database file from which each hit originates.
#'
#' @examples
#' hmm_file <- "your_hmm_file.hmm"
#' db_files <- c("database1.fasta", "database2.fasta")
#' search_results <- search_hmm(hmm_file, db_files, cpu = 2, incE = 1e-5)
#' @export
search_hmm <- function(hmm, dbs, cpu = 1, incE = 1e-6, hmmer_path = "") {
  tblout <- tempfile(fileext = ".tblout")
  if (is.null(names(dbs))) {
    names(dbs) <- dbs
  }
  hmmsearch_cmd <- if (nzchar(hmmer_path)) file.path(hmmer_path, "hmmsearch") else "hmmsearch"
  data.table::rbindlist(lapply(dbs, function(target) {
    system2("hmmsearch", c("--tblout", tblout, "--cpu", cpu, "--incE", incE, hmm, target))
    read_hmmer_tblout(tblout)
  }), id = "SeqFile")
}

#' Compare OTU hits depths between query and control
#'
#' @param ps A phyloseq object containing OTU (Operational Taxonomic Units) counts and (optionally) taxonomic information.
#' @param query_tblout Character; the HMMER output for the query genes.
#' @param control_tblout Character; the HMMER output for the control genes.
#' @param linker Function that links OTU identifiers between tblout tables and the phyloseq object. The function must take two arguments: SeqFile and Target.
#' @param taxrank Optional; character specifying the taxonomic rank at which to aggregate hits. Valid options depend on the taxonomic ranks present in the `ps` object. If not specified, no aggregation is performed.
#' @param phyloseq Logical; if TRUE, returns a phyloseq object, otherwise returns a data.table object. Defaults to TRUE.
#'
#' @return If `phyloseq` is TRUE (the default), this function returns a phyloseq object that contains the relative OTU depths of query vs control genes (aggregated by taxrank if specified). If `phyloseq` is FALSE, it returns a data.table object containing the calculated depths for queries and controls separately.
#'
#' @description `get_hits_depths` function calculates the relative OTU depths between the query and control genes, optionally aggregating hits at a specified taxonomic rank.
#'
#' @details The function prunes the phyloseq object to include only taxa that are present in both the query and control datasets based on the `linker` function. If a taxonomic rank is specified, the function aggregates hits at that rank and recalculates taxa names. It also warns about queries not found in the control. It performs calculates the depth ratios between query and control conditions.
#' Two linker functions are provided: mags_linker() and contigs_linker().
#'
#' @examples
#' ```
#' # Assume ps is a phyloseq object with taxa, otu_table etc., query_tblout and control_tblout are loaded
#' # Define a simple linker function
#' linker_function <- function(seq_file, target) {
#'   paste(seq_file, target, sep = "_")
#' }
#' # Call the get_hits_depths function without aggregation
#' get_hits_depths_result <- get_hits_depths(ps, query_tblout, control_tblout, linker_function)
#'
#' # Call the get_hits_depths function with aggregation at genus level
#' get_hits_depths_aggregated <- get_hits_depths(ps, query_tblout, control_tblout, linker_function, taxrank = "Genus")
#' ```
get_hits_depths <- function(ps, query_tblout, control_tblout, linker, taxrank = NULL, phyloseq = TRUE) {
  # stopifnot(length(query_tblout) == 1)
  # stopifnot(length(control_tblout) == 1)
  query_ps <- prune_taxa(unique(linker(query_tblout$SeqFile, query_tblout$Target)), ps)
  control_ps <- prune_taxa(unique(linker(control_tblout$SeqFile, control_tblout$Target)), ps)
  if (!is.null(taxrank)) {
    query_ps <- tax_glom(query_ps, taxrank)
    taxa_names(query_ps) <- apply(tax_table(query_ps), 1, function(x) paste(na.omit(x), collapse = ";"))
    control_ps <- tax_glom(control_ps, taxrank)
    taxa_names(control_ps) <- apply(tax_table(control_ps), 1, function(x) paste(na.omit(x), collapse = ";"))
    query_otu <- if (taxa_are_rows(ps)) {
      otu_table(query_ps)
    } else {
      t(otu_table(query_ps))
    }
    control_otu <- if (taxa_are_rows(ps)) {
      otu_table(control_ps)
    } else {
      t(otu_table(control_ps))
    }
    queries_not_controlled <- setdiff(rownames(query_otu), rownames(control_otu))
    if (length(queries_not_controlled)) {
      warning(queries_not_controlled, " do not have any control hits")
    }
    mer <- merge(
      as.data.table(query_otu, keep.rownames = "ID"),
      as.data.table(control_otu, keep.rownames = "ID"),
      by = "ID",
    )
    if (isFALSE(phyloseq)) {
      return(mer)
    }
    res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
    setnames(res, names(res), sub("\\.x$", "", names(res)))
    m <- as.matrix(res, rownames.value = mer$ID)
    # XXX: set NaNs (0/0) to zero
    m[is.nan(m)] <- 0
    otu_table(query_ps) <- otu_table(m, taxa_are_rows = T)
    return(query_ps)
  } else if (taxa_are_rows(ps)) {
    mer <- merge(
      as.data.table(t(colSums(otu_table(query_ps))))[, ID := "sp1"],
      as.data.table(t(colSums(otu_table(control_ps))))[, ID := "sp1"],
      by = "ID"
    )
    if (isFALSE(phyloseq)) {
      return(mer)
    }
    res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
    setnames(res, names(res), sub("\\.x$", "", names(res)))
    m <- as.matrix(res, rownames.value = mer$ID)
    # XXX: set NaNs (0/0) to zero
    m[is.nan(m)] <- 0
    return(phyloseq(otu_table(m, taxa_are_rows = T), access(ps, "sam_data"), access(ps, "phy_tree"), access(ps, "ref_seq")))
  } else {
    # taxa are cols
    mer <- merge(
      as.data.table(t(rowSums(otu_table(query_ps))))[, ID := "sp1"],
      as.data.table(t(rowSums(otu_table(control_ps))))[, ID := "sp1"],
      by = "ID"
    )
    if (isFALSE(phyloseq)) {
      return(mer)
    }
    res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
    setnames(res, names(res), sub("\\.x$", "", names(res)))
    m <- t(as.matrix(res, rownames.value = mer$ID))
    # XXX: set NaNs (0/0) to zero
    m[is.nan(m)] <- 0
    return(phyloseq(otu_table(m, taxa_are_rows = F), access(ps, "sam_data"), access(ps, "phy_tree"), access(ps, "ref_seq")))
  }
}

#' @export
mags_linker <- function(SeqName, Target) SeqName
#' @export
contigs_linker <- function(SeqName, Target) paste(sub("_[^_]+$", "", Target), SeqName, sep = "@")
