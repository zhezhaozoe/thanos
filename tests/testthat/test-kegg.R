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
