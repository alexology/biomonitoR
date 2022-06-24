test_that("bioco", {
  load(system.file("testdata", "ibmr_results.rda", package="biomonitoR"))
  load(system.file("testdata", "ibmr_taxa.rda", package="biomonitoR"))

  oglio_asb <- as_biomonitor(ibmr_taxa, group = "mf")
  oglio_agg <- aggregate_taxa(oglio_asb)

  to_remove <- c("Spirogyra dubia",
                 "Spirogyra gracilis",
                 "Spirogyra majuscola",
                 "Spirogyra parva")

  res <- ibmr(oglio_agg, coverage_coef = TRUE, exceptions = to_remove, traceB = TRUE)

  ibmr_expected <- ibmr_results$result
  names(ibmr_expected) <- paste("sito", 1:9, sep = "_")
  expect_equal(res$results, ibmr_expected)


  expect_message(ibmr(oglio_agg, coverage_coef = TRUE, exceptions = to_remove, traceB = TRUE), "Exceptions found, please see traceB for further information")


  # try with coverage classes
  ibmr_taxa_abu <- ibmr_taxa
  ibmr_taxa_abu[ibmr_taxa_abu == 1] <- 0.05
  ibmr_taxa_abu[ibmr_taxa_abu == 2] <- 0.5
  ibmr_taxa_abu[ibmr_taxa_abu == 3] <- 6
  ibmr_taxa_abu[ibmr_taxa_abu == 4] <- 25
  ibmr_taxa_abu[ibmr_taxa_abu == 5] <- 75

  oglio_asb_abu <- as_biomonitor(ibmr_taxa_abu, group = "mf")
  oglio_agg_abu <- aggregate_taxa(oglio_asb_abu)

  to_remove <- c("Spirogyra dubia",
                 "Spirogyra gracilis",
                 "Spirogyra majuscola",
                 "Spirogyra parva")

  res_abu <- ibmr(oglio_agg_abu, coverage_coef = FALSE, exceptions = to_remove, traceB = TRUE)

  expect_equal(res_abu$results, ibmr_expected)

  expect_equal(abundance_classes(ibmr_taxa_abu), ibmr_taxa)

})
