test_that("as_data_frame", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)

  macro_ex_wrong <- macro_ex
  macro_ex_wrong$Taxa <- as.character(macro_ex_wrong$Taxa)
  macro_ex_wrong[macro_ex_wrong$Taxa %in% "Acentrella", "Taxa"] <- "Acentrela"
  data_bio_wrong <- as_biomonitor(macro_ex_wrong, group = "mi", traceB = TRUE)
  out_put <- data.frame(excluded = "Acentrela", suggested = "Acentrella")

  expect_error(as.data.frame(data_bio, object = 2), "data_bio is of length 1")
  expect_equal(as.data.frame(data_bio_wrong, object = 2), out_put)
})


test_that("subset", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  subset_hydro <- subset(data_bio, taxa = "Hydropsychidae")
  subset_hydro_man <- as.data.frame(data_bio)
  subset_hydro_man <- subset_hydro_man[subset_hydro_man$Family %in% "Hydropsychidae", ]
  rownames(subset_hydro_man) <- NULL

  exclude_hydro <- subset(data_bio, exclude = "Hydropsychidae")
  exclude_hydro_man <- as.data.frame(data_bio)
  exclude_hydro_man <- exclude_hydro_man[!exclude_hydro_man$Family %in% "Hydropsychidae", ]
  rownames(exclude_hydro_man) <- NULL

  both_hydro <- subset(data_bio, taxa = "Hydropsychidae", exclude = "Cheumatopsyche")
  both_hydro_man <- as.data.frame(data_bio)
  both_hydro_man <- both_hydro_man[ both_hydro_man$Genus %in% "Hydropsyche", ]
  rownames(both_hydro_man) <- NULL

  both_hydro2 <- subset(data_bio, taxa = "Hydropsychidae", exclude = "Acentrella")

  expect_error(subset(data_bio, taxa = NULL, exclude = NULL), "taxa and exclude cannot be both null")
  expect_equal(as.data.frame(subset_hydro), subset_hydro_man)
  expect_equal(as.data.frame(exclude_hydro), exclude_hydro_man)
  expect_equal(as.data.frame(both_hydro), both_hydro_man)
  expect_equal(as.data.frame(both_hydro2), subset_hydro_man)
  expect_error(subset(data_bio, taxa = "Hydropsychide"),
               "None of the taxa provided were found in the data_bio database")
  expect_message(subset(data_bio, taxa = c("Hydropsychidae", "Hydropsychide")),
               "The following taxa were not find in the data_bio database and has been excluded: Hydropsychide")

  expect_message(subset(data_bio, exclude = c("Hydropsychidae", "Hydropsychide")),
                 "The following taxa were not find in the data_bio database and has been excluded: Hydropsychide")

  expect_message(subset(data_bio, taxa = c("Acentrella", "Acentrela"), exclude = c("Hydropsychidae", "Hydropsychide")),
                 "The following taxa were not find in the data_bio database and has been excluded: Acentrela, Hydropsychide")


})
