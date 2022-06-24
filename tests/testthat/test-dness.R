test_that("dness", {

  load(system.file("testdata", "tax_tree.rda", package="biomonitoR"))
  load(system.file("testdata", "comm_dness.rda", package="biomonitoR"))

  data(oglio)
  data_bio <- as_biomonitor(oglio, group = "mf", FUN = bin)
  data_agr <- aggregate_taxa(data_bio)
  data_dness <- dness(data_agr, method = "delta", complete = TRUE)
  data_dness_star <- dness(data_agr, method = "delta.st", complete = TRUE)
  data_dness_bin <- dness(data_agr, method = "delta.bin", complete = TRUE)

  taxdis <- vegan::taxa2dist(tax_tree, varstep=FALSE)
  veg_res <- vegan::taxondive(comm_dness, taxdis, match.force = FALSE)

  expect_identical(data_dness, veg_res$D, tolerance = 0.000001)
  expect_identical(data_dness_star, veg_res$Dstar, tolerance = 0.000001)
  expect_identical(data_dness_bin, veg_res$Dplus, tolerance = 0.000001)

})
