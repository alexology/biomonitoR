test_that("spain", {
  load(system.file("testdata", "aspt_bmwp_spain.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(aspt_bmwp_spain)
  data_agr <- aggregate_taxa(data_bio)
  spa_bmwp <- c(35, 18)
  spa_bmwp_ex <- c(29, 12)
  names(spa_bmwp) <- names(spa_bmwp_ex) <- c("Sample_1", "Sample_2")
  spa_aspt <- spa_bmwp / c(6, 3)
  tb_mes <- data.frame(Child = "Ferrissia", Parent = "Planorbidae")
  expect_equal(suppressMessages(bmwp(data_agr, method = "spa")),  spa_bmwp)
  expect_equal(suppressMessages(aspt(data_agr, method = "spa")),  spa_aspt)
  expect_message(bmwp(data_agr, method = "spa"), "Parent-child pairs found, please see traceB for further information" )
  expect_message(aspt(data_agr, method = "spa"), "Parent-child pairs found, please see traceB for further information" )
  expect_equal(suppressMessages(bmwp(data_agr, method = "spa", traceB = TRUE))$parent_child_pairs,  tb_mes)
  expect_equal(suppressMessages(aspt(data_agr, method = "spa", traceB = TRUE))$parent_child_pairs,  tb_mes)
  expect_equal(length(suppressMessages(bmwp(data_agr, method = "spa", traceB = TRUE))),  5)
  expect_equal(length(suppressMessages(aspt(data_agr, method = "spa", traceB = TRUE))),  5)
  expect_equal(suppressMessages(bmwp(data_agr, method = "spa", exceptions = "Ferrissia")),  spa_bmwp_ex)
  expect_equal(length(suppressMessages(bmwp(data_agr, method = "spa", exceptions = "Ferrissia", traceB = TRUE))),  5)
})


test_that("italy", {
  load(system.file("testdata", "aspt_bmwp_spain.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(aspt_bmwp_spain)
  data_agr <- aggregate_taxa(data_bio)
  ita_bmwp <- c(35, 21)
  names(ita_bmwp) <- c("Sample_1", "Sample_2")
  ita_aspt <- ita_bmwp / c(6, 4)
  ita_bmwp_agg <- c(21, 21)
  names(ita_bmwp_agg) <- c("Sample_1", "Sample_2")
  ita_aspt_agg <- ita_bmwp_agg / c(4, 4)
  expect_equal(suppressMessages(bmwp(data_agr, method = "ita")),  ita_bmwp)
  expect_equal(suppressMessages(aspt(data_agr, method = "ita")),  ita_aspt)
  expect_equal(suppressMessages(bmwp(data_agr, method = "ita", agg = TRUE)),  ita_bmwp_agg)
  expect_equal(suppressMessages(aspt(data_agr, method = "ita", agg = TRUE)),  ita_aspt_agg)
  expect_equal(length(suppressMessages(bmwp(data_agr, method = "ita", traceB = TRUE))),  5)
  expect_equal(length(suppressMessages(aspt(data_agr, method = "ita", traceB = TRUE))),  5)
  expect_equal(suppressMessages(bmwp(data_agr, method = "ita", agg = TRUE, traceB = TRUE))$composite_taxa,  c("Apataniidae", "Glossosomatidae"))
  expect_equal(suppressMessages(aspt(data_agr, method = "ita", agg = TRUE, traceB = TRUE))$composite_taxa,  c("Apataniidae", "Glossosomatidae"))

})


test_that("uk", {
  load(system.file("testdata", "aspt_bmwp_spain.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(aspt_bmwp_spain)
  data_agr <- aggregate_taxa(data_bio)
  uk_bmwp <- c(35, 21)
  names(uk_bmwp) <- c("Sample_1", "Sample_2")
  uk_aspt <- uk_bmwp / c(6, 4)
  uk_bmwp_agg <- c(21, 21)
  names(uk_bmwp_agg) <- c("Sample_1", "Sample_2")
  uk_aspt_agg <- uk_bmwp_agg / c(4, 4)
  expect_equal(suppressMessages(bmwp(data_agr, method = "uk")),  uk_bmwp)
  expect_equal(suppressMessages(aspt(data_agr, method = "uk")),  uk_aspt)
  expect_equal(suppressMessages(bmwp(data_agr, method = "uk", agg = TRUE)),  uk_bmwp_agg)
  expect_equal(suppressMessages(aspt(data_agr, method = "uk", agg = TRUE)),  uk_aspt_agg)
  expect_equal(length(suppressMessages(bmwp(data_agr, method = "uk", traceB = TRUE))),  5)
  expect_equal(length(suppressMessages(aspt(data_agr, method = "uk", traceB = TRUE))),  5)
  expect_equal(suppressMessages(bmwp(data_agr, method = "uk", agg = TRUE, traceB = TRUE))$composite_taxa,  c("Apataniidae", "Glossosomatidae"))
  expect_equal(suppressMessages(aspt(data_agr, method = "uk", agg = TRUE, traceB = TRUE))$composite_taxa,  c("Apataniidae", "Glossosomatidae"))
})




test_that("ita_bin", {
  load(system.file("testdata", "aspt_bmwp_spain.rda", package="biomonitoR"))
  data_bio <- suppressWarnings(as_biomonitor(aspt_bmwp_spain, FUN = bin))
  data_agr <- aggregate_taxa(data_bio)
  ita_bmwp <- c(35, 21)
  names(ita_bmwp) <- c("Sample_1", "Sample_2")
  ita_aspt <- ita_bmwp / c(6, 4)
  ita_bmwp_agg <- c(21, 21)
  names(ita_bmwp_agg) <- c("Sample_1", "Sample_2")
  ita_aspt_agg <- ita_bmwp_agg / c(4, 4)
  expect_equal(suppressMessages(bmwp(data_agr, method = "ita")),  ita_bmwp)
  expect_equal(suppressMessages(aspt(data_agr, method = "ita")),  ita_aspt)
})


test_that("ita_custom", {
  load(system.file("testdata", "aspt_bmwp_spain.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  ita_scores <- show_scores(index = "aspt", method = "ita")
  ita_scores_agg <- show_scores(index = "aspt", method = "uk")

  expect_equal(bmwp(data_agr, method = ita_scores_agg$scores, agg = ita_scores_agg$aggregation_rule),  bmwp(data_agr, method = "ita"))
  expect_equal(aspt(data_agr, method = ita_scores_agg$scores, agg = ita_scores_agg$aggregation_rule),  aspt(data_agr, method = "ita"))

})

