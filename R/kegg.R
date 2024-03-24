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
#' @importFrom msa msa
#' @export
get_kegg_msa <- function(ko, nmax = Inf, method = "Muscle", ...) {
  stopifnot(length(ko) == 1)
  genes <- names(KEGGREST::keggFind("genes", ko))
  if (nmax < Inf) {
    genes <- genes[1:min(length(genes), nmax)]
  }
  # keggGet expects at most 10 genes
  batches <- split(genes, 0:(length(genes)-1) %/% 10)
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
