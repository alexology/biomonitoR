test_that("bioco", {
  load(system.file("testdata", "data_uk_bioco.rda", package="biomonitoR"))
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(data_uk_bioco, dfref = ref_def)
  data_bio_bin <- suppressWarnings(as_biomonitor(data_uk_bioco, dfref = ref_def, FUN = bin))

  alien <- c("Potamopyrgus antipodarum", "Crangonyx pseudogracilis", "Crangonyx", "Dikerogammarus",
             "Dikerogammarus haemobaphes", "Chelicorophium curvispinum", "Corbicula fluminea", "Hypania invalida",
             "Dreissena polymorpha", "Dugesia tigrina")

  data_agr <- aggregate_taxa(data_bio)
  data_agr_bin <- aggregate_taxa(data_bio_bin)
  data_agr_alien <- aggregate_taxa(subset(data_bio, taxa = alien))

  data_rich <- suppressWarnings(richness(data_agr, tax_lev = "Family"))
  data_rich_alien <- richness(data_agr_alien, tax_lev = "Family")
  rich_ratio <- round(data_rich_alien/data_rich, 3)
  names(rich_ratio) <- NULL

  data_abu <- abundance(data_agr, tax_lev = "Taxa")
  data_abu_alien <- abundance(data_agr_alien, tax_lev = "Taxa")
  abu_ratio <- round(data_abu_alien/data_abu, 3)
  names(abu_ratio) <- NULL


  bioco_res <- bioco(data_agr, alien = alien, dfref = ref_def, digits = 3)
  bioco_res_zero <- bioco_res
  i <- unlist(lapply(bioco_res_zero, is.numeric))
  bioco_res_zero[i] <- 0
  bioco_res_bin <- bioco_res
  bioco_res_bin$aci <- bioco_res_bin$sbci <- NA


  expect_equal(bioco_res$rci, rich_ratio)
  expect_equal(bioco_res$aci, abu_ratio)

  expect_error(bioco(data_agr, alien = alien), "Please set dfref")
  expect_error(bioco(data_agr, dfref = "mi",alien = alien), "Please provide the dfref you used for as_biomonitor")
  expect_error(bioco(data_agr, dfref = ref_def), "Please provide a vector containing the names of alien taxa")
  expect_warning(bioco(data_agr, dfref = ref_def, alien = c("Potamopyrgus", "Ergo", "Agi")), "The following taxa are absent from the reference database: Ergo, Agi")
  expect_warning(bioco(data_agr, dfref = ref_def, alien = c("Ergo", "Agi")), "Ergo, Agi are not part of the taxonomic tree of data_agr")
  expect_equal(suppressWarnings(bioco(data_agr, dfref = ref_def, alien = c("Ergo", "Agi"))), bioco_res_zero)
  expect_equal(bioco(data_agr_bin, alien = alien, dfref = ref_def, digits = 3), bioco_res_bin)
})
