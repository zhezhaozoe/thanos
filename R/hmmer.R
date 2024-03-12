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
  msa(sequences, ...)
}

#' @importFrom Biostrings writeXStringSet unmasked
build_hmm <- function(aln) {
  faa <- tempfile(fileext = ".faa")
  writeXStringSet(unmasked(aln), file=faa)
  sto <- tempfile(fileext = ".sto")
  hmm <- tempfile(fileext = ".hmm")
  system2("esl-reformat", c("stockholm", faa), stdout = sto)
  system2("hmmbuild", c(hmm, sto))
  return(hmm)
}

#' @param ignore_target Should be TRUE for mag depths, FALSE for contig depths
#'
#' @import data.table
search_hmm <- function(hmm, dbs, ignore_target, cpu = 1, incE = 1e-6, target_sub_pattern = "_[^_]*$", target_sub_replacement = "") {
  tblout <- tempfile(fileext = ".tblout")
  hits <- rbindlist(lapply(dbs, function(target) {
    system2("hmmsearch", c("--tblout", tblout, "--cpu", cpu, "--incE", incE, hmm, target))
    read_hmmer_tblout(tblout)
  }), id = "SeqFile")
  if (ignore_target) {
    hits[, c("Target", "SeqFile") := .(SeqFile, NULL)]
  } else {
    hits[, c("Target", "SeqFile") := .(paste(sub(target_sub_pattern, target_sub_replacement, Target), SeqFile, sep = "@"), NULL)]
  }
}

#' @export
#' @import data.table
get_hits_depths <- function(ps, hits_hmm, controls_hmm, by_taxon = NULL, ...) {
  # aln <- do.call(get_kegg_msa, c(list(ko), msa_args))
  # hmm_ko <- build_hmm(aln)
  # hmm_control <- system.file("extdata", "DNGNGWU00001.hmm", package = "zzthanos")
  ignore_target <- is.null(attr(ps, "contigs"))
  hits_tblout <- search_hmm(hits_hmm, attr(ps, "seq_files"), ignore_target, ...)
  controls_tblout <- search_hmm(controls_hmm, attr(ps, "seq_files"), ignore_target, ...)
  if (!is.null(by_taxon)) {
    hits_depths <- as.data.table(otu_table(ps)[hits_tblout$Target, ])[, lapply(.SD, sum), by = list(Taxon = as.data.table(tax_table(ps)[hits_tblout$Target, ]@.Data)[[by_taxon]])]
    controls_depths <- as.data.table(otu_table(ps)[controls_tblout$Target, ])[, lapply(.SD, sum), by = list(Taxon = as.data.table(tax_table(ps)[controls_tblout$Target, ]@.Data)[[by_taxon]])]
    mer <- merge(hits_depths, controls_depths, by = "Taxon")
    res <- mer[, .SD, .SDcols = patterns("*\\.x$")] / mer[, .SD, .SDcols = patterns("*\\.y$")]
    setnames(res, names(res), sub("\\.x$", "", names(res)))
    melt(
      cbind(
        KO = ko,
        Taxon = mer$Taxon,
        res
      ),
      id.vars = c("KO", "Taxon"),
      variable.name = "Sample",
      value.name = "RelativeDepth"
    )
  } else {
    hits_depths <- colSums(otu_table(ps)[hits_tblout$Target, ])
    controls_depths <- colSums(otu_table(ps)[controls_tblout$Target, ])
    res <- hits_depths / controls_depths
    data.table(
      KO = ko,
      Taxon = NA,
      Sample = names(res),
      RelativeDepth = res
    )
  }
}
