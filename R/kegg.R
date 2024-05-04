#' Retrieve and Align KEGG Gene Sequences
#'
#' This function retrieves amino acid sequences for genes associated with a given KEGG Orthology (KO) identifier and performs a multiple sequence alignment (MSA) on them.
#'
#' @param ko A character string specifying the KEGG Orthology (KO) identifier.
#' @param nmax The maximum number of genes to include in the alignment (Default: no limit).
#' @param ... Additional arguments passed to the `msa` function for performing the multiple sequence alignment.
#'
#' @return An object of class `msa` representing the multiple sequence alignment of the retrieved sequences. The name of the KO identifier (or, if NULL, the KO identifier itself) is stored as an attribute named \code{"name"} of the returned object.
#'
#' @examples
#' \dontrun{
#'   aln <- get_kegg_msa("K02533")
#'   plot(aln)
#' }
#'
#' @importFrom KEGGREST keggFind keggGet
#' @import msa
#' @export
get_kegg_msa <- function(ko, nmax = Inf, method = "Muscle", ...) {
  stopifnot(length(ko) == 1)
  genes <- names(KEGGREST::keggFind("genes", ko))
  if (nmax < Inf) {
    genes <- genes[seq_len(min(length(genes), nmax))]
  }
  # keggGet expects at most 10 genes
  batches <- split(genes, 0:(length(genes) - 1) %/% 10)
  sequence_sets <- sapply(batches, function(batch) {
    Sys.sleep(.1)
    KEGGREST::keggGet(batch, "aaseq")
  })
  sequences <- Reduce(`c`, sequence_sets)
  aln <- msa::msa(sequences, method = method, ...)
  attr(aln, "name") <- if (!is.null(names(ko))) names(ko) else ko
  aln
}

#' Retrieve KEGG Orthology IDs (KOs) associated with a KEGG module
#'
#' This function fetches KEGG Orthology IDs (KOs) that are associated with a
#' specified KEGG module. It queries the KEGG database for the module, extracts
#' the KOs, and returns them as a vector.
#'
#' @param module A character string specifying the KEGG module ID (e.g., "M00001").
#'               Only one module ID should be provided.
#'
#' @return A character vector of KEGG Orthology IDs (KOs) associated with the
#'         specified module. If multiple KOs are associated, they are returned
#'         as separate elements in the vector.
#'
#' @examples
#' kegg_kos <- get_kegg_kos_from_module("M00001")
#' print(kegg_kos)
#'
#' @export
#'
#' @importFrom KEGGREST keggGet
#'
#' @note: The function relies on the `KEGGREST` package for fetching data
#' from the KEGG database. Ensure that this package is installed.
#' into your R session.
get_kegg_kos_from_module <- function(module) {
  stopifnot(length(module) == 1)
  m <- KEGGREST::keggGet(module)[[1]]
  kos_raw <- names(m$ORTHOLOGY)
  kos <- unique(unlist(strsplit(kos_raw, "[,+]")))
  # keggGet expects at most 10 genes
  batches <- split(kos, 0:(length(kos) - 1) %/% 10)
  symbols_list <- sapply(batches, function(batch) {
    Sys.sleep(.1)
    sapply(KEGGREST::keggGet(batch), `[[`, "SYMBOL")
  })
  setNames(kos, Reduce(`c`, symbols_list))
}

#' Retrieve KEGG reactions from a specified module as a graph
#'
#' This function fetches reactions for a specified KEGG module using the
#' KEGGREST API and processes them to create a graph representation. In
#' this graph, the nodes are compounds and the edges are reactions.

#' @param module A character string specifying the KEGG module ID (e.g.
#' "M00001").
#'
#' @details The KEGG graph has compounds as nodes and reactions as
#' edges. Each edge can actually be associated with multiple reactions,
#' in which case all of them occur. (Typically, they are the same
#' chemical transformation with different cofactors or performed by
#' different enzymes.) In turn, each reaction can be associated with one
#' or more KOs. The same KO can be associated to multiple reactions as
#' well.
#'
#' @return A data.table object representing the simplified graph with
#' the following columns: `from` (ID of the starting compound), `to`
#' (ID of the resulting compound), `reaction` (concatenated string of
#' reaction IDs and their associated KOs, separated by "|"), `from_name`
#' (name of the starting compound), and `to_name` (name of the resulting
#' compound).
#'
#' @examples
#' \dontrun{
#'   module <- "M00001"
#'   graph <- get_kegg_reactions_from_module(module)
#'   print(graph)
#' }
#'
#' @import data.table
#' @importFrom KEGGREST keggGet
#' @export
get_kegg_reactions_from_module <- function(module) {
  stopifnot(length(module) == 1)
  m <- KEGGREST::keggGet(module)[[1]]
  graph <- setNames(setDT(tstrsplit(m$REACTION, " -> ")), c("from", "to"))
  orthology <- gsub(".*\\[RN:", "", m$ORTHOLOGY)
  # Stop if there is a comma: our assumption may be wrong
  stopifnot(!any(grepl(",", orthology)))
  orthology <- strsplit(gsub("\\]", "", orthology), " ")
  inverted_orthology <- inverted_names(orthology)
  if (is.list(inverted_orthology)) {
    inverted_orthology <- sapply(inverted_orthology, paste, collapse = ",")
  }
  reactions <- sapply(strsplit(names(m$REACTION), ","), function(reactions) {
    paste(reactions, inverted_orthology[reactions], sep = ":", collapse = "|")
  })
  graph$reaction <- reactions
  graph <- graph[, .(to = strsplit(to, " \\+ ")[[1]]), by = c("from", "reaction")]
  graph <- graph[, .(from = strsplit(from, " \\+ ")[[1]]), by = c("to", "reaction")]
  graph <- unique(graph)
  graph$from_name <- m$COMPOUND[graph$from]
  graph$to_name <- m$COMPOUND[graph$to]
  graph
}
