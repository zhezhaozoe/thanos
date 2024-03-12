mag_depths_file <- "/g/scb/bork/marotta/thanos/GenomeBinning/depths/bins/bin_depths_summary.tsv"
seq_files <- Sys.glob("/g/scb/bork/marotta/thanos/Annotation/Prokka/*.faa")
names(seq_files) <- sub("\\.faa", "", basename(seq_files))
samplesheet_file <- "/g/bork3/home/marotta/samplesheet_v2.csv"
gtdb_file <- "/g/scb/bork/marotta/thanos/Taxonomy/GTDB-Tk/gtdbtk_summary.tsv"

test_that("we can import mag depths file", {
  samplesheet <- sample_data(read.csv(samplesheet_file, header = T, row.names = 1))
  gtdb <- tax_table(read_gtdbtk(gtdb_file))
  ps <- import_mag_depths(mag_depths_file, samplesheet, gtdb, seq_files)
})

test_that("we can run hmmer", {
  ps <- import_mag_depths(mag_depths_file, samplesheet, gtdb, seq_files)
  aln <- get_kegg_msa("K02588", method = "Muscle")
  hmm_ko <- build_hmm(aln)
  hmm_control <- system.file("extdata", "DNGNGWU00001.hmm", package = "zzthanos")
  search_hmm(hmm_ko, attr(ps, "seq_files"), ignore_target = T)
})
