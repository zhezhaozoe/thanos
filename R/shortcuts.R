tblout_from_ko <- function(ko, method = "Muscle", cpu = 1, incE = 1e-6) {
  aln <- get_kegg_msa(ko, method = method)
  hmm <- build_hmm(aln)
  search_hmm(hmm, dbs, cpu = cpu, incE = incE)
}

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

get_hits_depths_from_kos <- function(kos, ps, dbs, control_tblout, linker, msa_method = "Muscle", cpu = 1, incE = 1e-6, taxrank = NULL) {
  lapply(kos, function(ko) {
    query_tblout <- tblout_from_ko(ko, method = msa_method, cpu = cpu, incE = incE)
    get_hits_depths(ps, query_tblout, control_tblout, linker, taxrank = taxrank)
  })
}

get_contigs_hits_depths_from_kos <- function(kos, depths_files, pattern, replacement, dbs, control_tblout, linker = contigs_linker, msa_method = "Muscle", cpu = 1, incE = 1e-6, verbose = F) {
  lapply(kos, function(ko) {
    query_tblout <- tblout_from_ko(ko, method = msa_method, cpu = cpu, incE = incE)
    get_contigs_hits_depths(depths_files, pattern, replacement, query_tblout, control_tblout, linker = linker, verbose = verbose)
  })
}
