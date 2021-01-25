test_that("life", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))
  load(system.file("testdata", "metrics_uk_data.rda", package="biomonitoR"))
  load(system.file("testdata", "taxa_uk_data.rda", package="biomonitoR"))
  load(system.file("testdata", "to_change.rda", package="biomonitoR"))

  taxa_uk_data.asb <- suppressMessages(as_biomonitor(taxa_uk_data, dfref = ref_def, to_change = to_change, traceB = FALSE))
  taxa_uk_data.agg <- aggregate_taxa(taxa_uk_data.asb)
  res_life <- round(life(taxa_uk_data.agg, method = "life_2017"), 2)
  res_whpt <- round(whpt(taxa_uk_data.agg, method = "uk", agg = TRUE), 2)
  res_psi <- round(psi(taxa_uk_data.agg, method = "extence"), 2)
  res_whpt_bmwp <- round(whpt(taxa_uk_data.agg, method = "uk", metric = "bmwp", agg = TRUE), 2)
  names(res_life) <- names(res_whpt) <- names(res_whpt_bmwp) <- names(res_psi) <- NULL
  expect_equal(res_life, metrics_uk_data$LIFE_FAMILY_INDEX)
  expect_equal(res_whpt, metrics_uk_data$WHPT_ASPT)
  expect_equal(res_whpt_bmwp, metrics_uk_data$WHPT_TOTAL)
  expect_equal(res_psi, metrics_uk_data$PSI_FAMILY_SCORE)
})
