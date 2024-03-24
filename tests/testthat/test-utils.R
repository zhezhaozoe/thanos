test_that("Shortening the species names works", {
  expect_equal(shorten_species_name("Canis lupus"), "C. lupus")
  expect_equal(shorten_species_name("Lynx"), "Lynx")
  expect_equal(shorten_species_name("Candidatus Mycoplasma hyopneumoniae"), "Ca. M. hyopneumoniae")
  expect_equal(shorten_species_name("Proteobacteria species incertae sedis"), "Proteobacteria sp. inc. s.")
})
