test_that("bio_indices", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)

  #ept
  expect_equal(ept(data_agr), get_taxa_richness(data_agr, taxa = c("Ephemeroptera", "Plecoptera", "Trichoptera"), tax_lev = "Taxa"))
  expect_equal(ept(data_agr, tax_lev = "Family"), get_taxa_richness(data_agr, taxa = c("Ephemeroptera", "Plecoptera", "Trichoptera"), tax_lev = "Family"))

  #igold
  expect_equal(igold(data_agr), 1 - get_taxa_abundance(data_agr, taxa = c("Oligochaeta", "Diptera", "Gastropoda"), rel = TRUE))

  # eptd
  taxa_eptd <- c("Heptageniidae", "Leptophlebiidae",  "Polycentropodidae", "Empididae")
  expect_equal(eptd(data_agr), log(get_taxa_abundance(data_agr, taxa = taxa_eptd) + 1, 10))
  expect_equal(eptd(data_agr, base = 2), log(get_taxa_abundance(data_agr, taxa = taxa_eptd) + 1, 2))
  expect_equal(eptd(data_agr, eptd_families = taxa_eptd, base = 2), log(get_taxa_abundance(data_agr, taxa = taxa_eptd) + 1, 2))
})


