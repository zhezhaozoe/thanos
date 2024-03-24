#' Test for get_kegg_kos_from_module
#'
#' Tests the get_kegg_kos_from_module function to ensure it properly retrieves KEGG KOs from the specified module.
#'
#' @test
test_that("get_kegg_kos_from_module retrieves KOs correctly", {
  results <- list(
    M00099 = c("K00654", "K04708", "K04709", "K04710", "K24621", "K24622", "K23727", "K04712", "K01441", "K12348", "K12349")
  )
  lapply(seq_along(results), function(i) {
    expect_equal(
      get_kegg_kos_from_module(names(results)[i]),
      results[[i]]
    )
  })
})

test_that("get_kegg_msa stops for vector length not equal to 1", {
  expect_error(get_kegg_msa(c("K00010", "K00020")))
})

test_that("get_kegg_msa returns the correct object", {
  ko <- "K00010"
  nmax <- 10
  result <- get_kegg_msa(ko, nmax = nmax)
  expect_equal(dim(result)[1], nmax)
  expect_s4_class(result, "MsaAAMultipleAlignment")
})

test_that("get_kegg_msa returns an object with the correct name attribute", {
  ko <- "K00010"
  result <- get_kegg_msa(ko, nmax = 10)
  expect_equal(attr(result, "name"), ko)

  ko <- c(geneA = "K00010")
  result <- get_kegg_msa(ko, nmax = 10)
  expect_equal(attr(result, "name"), "geneA")
})

