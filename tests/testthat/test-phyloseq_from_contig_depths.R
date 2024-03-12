contig_depths_files <- Sys.glob("/g/scb/bork/marotta/thanos/GenomeBinning/depths/contigs/*.txt.gz")
names(contig_depths_files) <- sub("MEGAHIT-(.*)-depth.txt.gz", "\\1", basename(contig_depths_files))
seq_files <- Sys.glob("/g/scb/bork/marotta/thanos/Annotation/Prodigal/MEGAHIT/group-*/*.faa.gz")
names(seq_files) <- sub("\\.faa.gz", "", basename(seq_files))
samplesheet_file <- "/g/bork3/home/marotta/samplesheet_v2.csv"

test_that("we can import depths from contigs", {
  samplesheet <- sample_data(read.csv(samplesheet_file, header = T, row.names = 1))
  ps <- import_contig_depths(contig_depths_files, samplesheet, sub_pattern = "MEGAHIT-group-\\d+-(.*).bam", sub_replacement = "\\1", seq_files)
})

test_that("we can run hmmer", {
  ps <- import_contig_depths(contig_depths_files, samplesheet, sub_pattern = "MEGAHIT-group-\\d+-(.*).bam", sub_replacement = "\\1", seq_files)
  aln <- get_kegg_msa("K02588", method = "Muscle")
  hmm_ko <- build_hmm(aln)
  hmm_control <- system.file("extdata", "DNGNGWU00001.hmm", package = "zzthanos")
  hits_ko <- search_hmm(hmm_ko, attr(ps, "seq_files"), ignore_target = F)
  hits_control <- search_hmm(hmm_control, attr(ps, "seq_files"), ignore_target = F)
})
