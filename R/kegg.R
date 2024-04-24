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
    genes <- genes[1:min(length(genes), nmax)]
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
  m <- KEGGREST::keggGet(module)
  kos_raw <- names(m[[1]]$ORTHOLOGY)
  unlist(strsplit(kos_raw, ","))
}

#' @import data.table
#' @importFrom KEGGREST keggGet
#' @export
get_kegg_reactions_from_module <- function(module) {
  # NOTE: make a simplification: for reactions with multiple input/output compounds, the multiple compounds are treated as just one. Let's see how the graph turns out before making anything more complex.
  # NOTE: in general, there is a many-to-many relationship between reaction IDs and KOs. Here we assume that the relationship is just one reaction ID -> many KOs.
  # NOTE: I now think the correct interface is to make a bipartite graph with reactions and compounds as nodes. the links can be only from compound to reaction and from reaction to another compound. there can be multiple reactions for the same link, in which case they can all occur, typically with different cofactors (e.g. with either ATP or ADP as the P donor). Then, independently of this, each KO is associated to one or more reactions.
  stopifnot(length(module) == 1)
  m <- KEGGREST::keggGet(module)[[1]]
  graph <- setNames(setDT(tstrsplit(m$REACTION, " -> ")), c("from", "to"))
  orthology <- gsub(".*\\[RN:", "", m$ORTHOLOGY)
  # Stop if there is a comma: our assumption may be wrong
  stopifnot(!any(grepl(",", orthology)))
  orthology <- strsplit(gsub("\\]", "", orthology), " ")
  inverted_orthology <- inverted_names(unlist(orthology))
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
