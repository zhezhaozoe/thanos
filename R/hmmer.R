#' @importFrom KEGGREST keggFind keggGet
#' @import msa
get_kegg_msa <- function(ko, ...) {
  genes <- names(keggFind("genes", ko))
  batches <- split(genes, 0:(length(genes)-1) %/% 10)
  sequence_sets <- sapply(batches[1:2], function(batch) {
    Sys.sleep(.1)
    keggGet(batch, "aaseq")
  })
  sequences <- Reduce(`c`, sequence_sets)
  aln <- msa(sequences, ...)
  attr(aln, "name") <- ko
  aln
}

#' @importFrom Biostrings writeXStringSet unmasked
build_hmm <- function(aln) {
  faa <- tempfile(fileext = ".faa")
  writeXStringSet(unmasked(aln), file=faa)
  sto <- tempfile(fileext = ".sto")
  hmm <- tempfile(fileext = ".hmm")
  extargs <- if (! is.null(attr(aln, "name"))) {
    c("-n", attr(aln, "name"))
  } else {
    c()
  }
  system2("esl-reformat", c("stockholm", faa), stdout = sto)
  system2("hmmbuild", c(extargs, hmm, sto))
  return(hmm)
}

#' @import data.table
#' @export
search_hmm <- function(hmm, dbs, cpu = 1, incE = 1e-6) {
  tblout <- tempfile(fileext = ".tblout")
  rbindlist(lapply(dbs, function(target) {
    system2("hmmsearch", c("--tblout", tblout, "--cpu", cpu, "--incE", incE, hmm, target))
    read_hmmer_tblout(tblout)
  }), id = "SeqFile")
}

#' @import data.table
#' @export
get_hits_depths <- function(ps, query_tblout, control_tblout, linker, taxrank = NULL) {
  query_ps <- prune_taxa(unique(linker(query_tblout$SeqFile, query_tblout$Target)), ps)
  control_ps <- prune_taxa(unique(linker(control_tblout$SeqFile, control_tblout$Target)), ps)
  if (! is.null(taxrank)) {
    query_ps <- tax_glom(query_ps, taxrank)
    taxa_names(query_ps) <- apply(tax_table(query_ps), 1, paste, collapse = ";_;")
    control_ps <- tax_glom(control_ps, taxrank)
    taxa_names(control_ps) <- apply(tax_table(control_ps), 1, paste, collapse = ";_;")
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
    mer <- merge(
      as.data.table(query_otu, keep.rownames = "ID"),
      as.data.table(control_otu, keep.rownames = "ID"),
      by = "ID",
      all.x = T
    )
    res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
    setnames(res, names(res), sub("\\.x$", "", names(res)))
    m <- as.matrix(res, rownames.value = mer$ID)
    otu_table(query_ps) <- otu_table(m, taxa_are_rows = T)
    return(query_ps)
  } else if (taxa_are_rows(ps)) {
    m <- matrix(
      colSums(otu_table(query_ps)) / colSums(otu_table(control_ps)),
      nrow = 1)
    colnames(m) <- colnames(otu_table(ps))
    return(phyloseq(otu_table(m, taxa_are_rows(ps)), access(ps, "sam_data"), access(ps, "phy_tree"), access(ps, "ref_seq")))
  }
}

#' @export
mags_linker <- function(SeqName, Target) SeqName
#' @export
contigs_linker <- function(SeqName, Target) paste(SeqName, sub("_[^_]+$", "", Target), sep = "@")
