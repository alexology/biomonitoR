test_that("scores_aspt", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))
  spa_scores <- show_scores(index = "aspt", method = "spa")$scores
  names(spa_scores)[1] <- "Taxa"
  spa_bio <- as_biomonitor(spa_scores, dfref = ref_def)
  expect_equal(length(spa_bio),  1)

  ita_scores <- show_scores(index = "aspt", method = "ita")
  ita_scores_1 <- ita_scores$scores
  ita_scores_2 <- data.frame(ita_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(ita_scores$aggregation_rule)))
  ita_scores_3 <- data.frame(ita_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(ita_scores$aggregation_rule)))
  names(ita_scores_1)[1] <- names(ita_scores_2)[1] <- names(ita_scores_3)[1] <- "Taxa"

  ita_bio_1 <- as_biomonitor(ita_scores_1, dfref = ref_def)
  ita_bio_2 <- as_biomonitor(ita_scores_2, dfref = ref_def)
  ita_bio_3 <- as_biomonitor(ita_scores_3, dfref = ref_def)

  expect_equal(length(ita_bio_1), 1)
  expect_equal(length(ita_bio_2), 1)
  expect_equal(length(ita_bio_3), 1)

  uk_scores <- show_scores(index = "aspt", method = "uk")
  uk_scores_1 <- uk_scores$scores
  uk_scores_2 <- data.frame(uk_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  uk_scores_3 <- data.frame(uk_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  names(uk_scores_1)[1] <- names(uk_scores_2)[1] <- names(uk_scores_3)[1] <- "Taxa"

  uk_bio_1 <- as_biomonitor(uk_scores_1, dfref = ref_def)
  uk_bio_2 <- as_biomonitor(uk_scores_2, dfref = ref_def)
  uk_bio_3 <- as_biomonitor(uk_scores_3, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)
  expect_equal(length(uk_bio_2), 1)
  expect_equal(length(uk_bio_3), 1)

})


test_that("scores_bmwp", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))
  spa_scores <- show_scores(index = "bmwp", method = "spa")$scores
  names(spa_scores)[1] <- "Taxa"
  spa_bio <- as_biomonitor(spa_scores, dfref = ref_def)
  expect_equal(length(spa_bio),  1)

  ita_scores <- show_scores(index = "bmwp", method = "ita")
  ita_scores_1 <- ita_scores$scores
  ita_scores_2 <- data.frame(ita_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(ita_scores$aggregation_rule)))
  ita_scores_3 <- data.frame(ita_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(ita_scores$aggregation_rule)))
  names(ita_scores_1)[1] <- names(ita_scores_2)[1] <- names(ita_scores_3)[1] <- "Taxa"

  ita_bio_1 <- as_biomonitor(ita_scores_1, dfref = ref_def)
  ita_bio_2 <- as_biomonitor(ita_scores_2, dfref = ref_def)
  ita_bio_3 <- as_biomonitor(ita_scores_3, dfref = ref_def)

  expect_equal(length(ita_bio_1), 1)
  expect_equal(length(ita_bio_2), 1)
  expect_equal(length(ita_bio_3), 1)

  uk_scores <- show_scores(index = "bmwp", method = "uk")
  uk_scores_1 <- uk_scores$scores
  uk_scores_2 <- data.frame(uk_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  uk_scores_3 <- data.frame(uk_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  names(uk_scores_1)[1] <- names(uk_scores_2)[1] <- names(uk_scores_3)[1] <- "Taxa"

  uk_bio_1 <- as_biomonitor(uk_scores_1, dfref = ref_def)
  uk_bio_2 <- as_biomonitor(uk_scores_2, dfref = ref_def)
  uk_bio_3 <- as_biomonitor(uk_scores_3, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)
  expect_equal(length(uk_bio_2), 1)
  expect_equal(length(uk_bio_3), 1)

})

test_that("scores_psi", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))

  uk_scores <- show_scores(index = "psi", method = "extence")
  uk_scores_1 <- uk_scores$fssr_groups
  names(uk_scores_1)[1] <- "Taxa"

  uk_bio_1 <- as_biomonitor(uk_scores_1, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)


})


test_that("scores_epsi", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))

  uk_scores <- show_scores(index = "epsi", method = "uk")
  uk_scores_1 <- uk_scores$scores
  uk_scores_2 <- data.frame(uk_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  uk_scores_3 <- data.frame(uk_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  names(uk_scores_1)[1] <- names(uk_scores_2)[1] <- names(uk_scores_3)[1] <- "Taxa"


  uk_bio_1 <- suppressWarnings(as_biomonitor(uk_scores_1, dfref = ref_def))
  uk_bio_2 <- as_biomonitor(uk_scores_2, dfref = ref_def)
  uk_bio_3 <- as_biomonitor(uk_scores_3, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)
  expect_equal(length(uk_bio_2), 1)
  expect_equal(length(uk_bio_3), 1)

})


test_that("scores_life_extence", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))

  uk_scores <- show_scores(index = "life", method = "extence")
  uk_scores_1 <- uk_scores$fs_groups
  uk_scores_2 <- data.frame(uk_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  uk_scores_3 <- data.frame(uk_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  names(uk_scores_1)[1] <- names(uk_scores_2)[1] <- names(uk_scores_3)[1] <- "Taxa"

  uk_bio_1 <- suppressWarnings(as_biomonitor(uk_scores_1, dfref = ref_def))
  uk_bio_2 <- as_biomonitor(uk_scores_2, dfref = ref_def)
  uk_bio_3 <- as_biomonitor(uk_scores_3, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)
  expect_equal(length(uk_bio_2), 1)
  expect_equal(length(uk_bio_3), 1)

})



test_that("scores_life_2017", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))

  uk_scores <- show_scores(index = "life", method = "life_2017")
  uk_scores_1 <- uk_scores$fs_groups
  uk_scores_2 <- data.frame(uk_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  uk_scores_3 <- data.frame(uk_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  names(uk_scores_1)[1] <- names(uk_scores_2)[1] <- names(uk_scores_3)[1] <- "Taxa"

  uk_bio_1 <- suppressWarnings(as_biomonitor(uk_scores_1, dfref = ref_def))
  uk_bio_2 <- as_biomonitor(uk_scores_2, dfref = ref_def)
  uk_bio_3 <- as_biomonitor(uk_scores_3, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)
  expect_equal(length(uk_bio_2), 1)
  expect_equal(length(uk_bio_3), 1)

})


test_that("scores_whpt", {
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))

  uk_scores <- show_scores(index = "whpt", method = "uk")
  uk_scores_1 <- uk_scores$scores
  uk_scores_2 <- data.frame(uk_scores$aggregation_rule[ , "Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  uk_scores_3 <- data.frame(uk_scores$aggregation_rule[ , "Correct_Taxon"], fake = rep(2, nrow(uk_scores$aggregation_rule)))
  names(uk_scores_1)[1] <- names(uk_scores_2)[1] <- names(uk_scores_3)[1] <- "Taxa"

  uk_bio_1 <- suppressWarnings(as_biomonitor(uk_scores_1, dfref = ref_def))
  uk_bio_2 <- as_biomonitor(uk_scores_2, dfref = ref_def)
  uk_bio_3 <- as_biomonitor(uk_scores_3, dfref = ref_def)

  expect_equal(length(uk_bio_1), 1)
  expect_equal(length(uk_bio_2), 1)
  expect_equal(length(uk_bio_3), 1)

})

