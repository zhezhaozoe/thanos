test_that("build_hmm handles different inputs correctly", {
  # Test for character input that does not exist
  expect_error(build_hmm("nonexistentfile.afa"), "File doesn't exist. Please provide a valid path")

  # Setup for MsaAAMultipleAlignment class input
  aln <- get_kegg_msa("K02588", nmax = 9)
  # Test for MsaAAMultipleAlignment input
  expect_silent({
    tmp <- build_hmm(aln)
    unlink(tmp)
  })

  # Test for unrecognized format input
  expect_error(build_hmm(list(a=1)), "Unrecognised format. aln can be either the path to an alignment in fasta format or an object of class from the msa package")
  
  # Test with afa path
  control_afa <- system.file("extdata", "controls", "bac120_r214_reps_PF01025.20.afa", package = "zzthanos")
  expect_silent({
    tmp <- build_hmm(control_afa)
    unlink(tmp)
  })
})
