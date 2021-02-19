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



test_that("remove_taxa", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)

  remove_hydroptilidae <- macro_ex[! macro_ex[, "Taxa"] %in% c("Hydroptilidae", "Hydroptila"), ]
  remove_hydroptila <- macro_ex[! macro_ex[, "Taxa"] %in% "Hydroptila", ]
  remove_hydroptila_dae <- macro_ex[! macro_ex[, "Taxa"] %in% c("Hydroptilidae", "Hydroptila"), ]
  remove_hydroptilidae_acentrella <- macro_ex[! macro_ex[, "Taxa"] %in% c("Hydroptilidae", "Hydroptila",  "Acentrella"), ]
  remove_hydrophilidae <- macro_ex[! macro_ex[, "Taxa"] %in% "Laccobius", ]

  remove_hydroptilidae$Taxa <- as.character(remove_hydroptilidae$Taxa)
  remove_hydroptila$Taxa <- as.character(remove_hydroptila$Taxa)
  remove_hydroptila_dae$Taxa <- as.character(remove_hydroptila_dae$Taxa)
  remove_hydroptilidae_acentrella$Taxa <- as.character(remove_hydroptilidae_acentrella$Taxa)
  remove_hydrophilidae$Taxa <- as.character(remove_hydrophilidae$Taxa)

  rownames(remove_hydroptilidae) <- rownames(remove_hydroptila) <- rownames(remove_hydroptila_dae) <-
    rownames(remove_hydroptilidae_acentrella) <- rownames(remove_hydrophilidae) <- NULL

  expect_equal(remove_taxa(data_agr, taxa = "Hydroptilidae"), remove_hydroptilidae)
  expect_equal(remove_taxa(data_agr, taxa = "Hydroptila"), remove_hydroptila)
  expect_equal(remove_taxa(data_agr, taxa = c("Hydroptilidae", "Hydroptila")), remove_hydroptila_dae)
  expect_equal(remove_taxa(data_agr, taxa = c("Hydroptilidae",  "Acentrella")), remove_hydroptilidae_acentrella)
  expect_equal(remove_taxa(data_agr, taxa = "Hydrophilidae"), remove_hydrophilidae)
  expect_error(remove_taxa(data_agr, taxa = NULL), "Please provide at least taxon name")
  expect_error(remove_taxa(data_agr, taxa = "ignita"), "None of the taxa provided were found in data_agr")
  expect_message(remove_taxa(data_agr, taxa = c("Acentrella","ignita", "ergo")), "The following taxa were not find in data_agr and has been excluded: ignita, ergo")
})
