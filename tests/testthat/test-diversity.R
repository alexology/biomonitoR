test_that("vegan", {
  data( macro_ex )
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_veg_taxa <- convert_to_vegan(data_agr, tax_lev = "Taxa")
  data_veg_gen <- convert_to_vegan(data_agr, tax_lev = "Genus")
  data_veg_fam <- convert_to_vegan(data_agr, tax_lev = "Family")
  expect_equal(shannon(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa))
  expect_equal(shannon(data_agr, tax_lev = "Genus"),  vegan::diversity(data_veg_gen))
  expect_equal(shannon(data_agr, tax_lev = "Taxa", base = 2 ),  vegan::diversity(data_veg_taxa, base = 2))
  expect_equal(simpson(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa, index = "simpson"))
  expect_equal(simpson(data_agr, tax_lev = "Genus"),  vegan::diversity(data_veg_gen, index = "simpson"))
  expect_equal(invsimpson(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa, index = "invsimpson"))
  expect_equal(invsimpson(data_agr, tax_lev = "Genus"),  vegan::diversity(data_veg_gen, index = "invsimpson"))
  expect_equal(shannon(data_agr, tax_lev = "Family"),  vegan::diversity(data_veg_fam))
  expect_equal(shannon(data_agr, tax_lev = "Genus"),  vegan::diversity(data_veg_gen))
  expect_equal(fisher(data_agr, tax_lev = "Taxa"),  vegan::fisher.alpha(data_veg_taxa), tolerance=1e-4)
  expect_equal(esimpson(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa, index = "invsimpson")/vegan::specnumber(data_veg_taxa))
})


test_that("vegan_bin", {
  data(macro_ex)
  data_bio <- suppressWarnings(as_biomonitor(macro_ex, FUN = bin))
  data_agr <- aggregate_taxa(data_bio)
  data_veg_taxa <- convert_to_vegan(data_agr, tax_lev = "Taxa")
  data_veg_fam <- convert_to_vegan(data_agr, tax_lev = "Family")
  expect_equal(shannon(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa))
  expect_equal(shannon(data_agr, tax_lev = "Taxa", base = 2 ),  vegan::diversity(data_veg_taxa, base = 2))
  expect_equal(simpson(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa, index = "simpson"))
  expect_equal(invsimpson(data_agr, tax_lev = "Taxa"),  vegan::diversity(data_veg_taxa, index = "invsimpson"))
  expect_equal(shannon(data_agr, tax_lev = "Family"),  vegan::diversity(data_veg_fam))
  #expect_equal(fisher(data_agr, tax_lev = "Taxa"),  vegan::fisher.alpha(data_veg_taxa), tolerance=1e-4)
})


test_that("past", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  # results compared with those of the past software version 4.04 accessed on 2021-01-11
  load(system.file("testdata", "past_results.rda", package="biomonitoR"))
  past_results <- as.data.frame (past_results[,-1], row.names = past_results[,1])
  past_mar <- past_results[, "Margalef"]
  past_brill <- past_results[, "Brillouin"]
  past_men <- past_results[, "Menhinick"]
  past_berpar <- past_results[, "Berger-Parker"]
  names(past_mar) <- names(past_brill) <- names(past_men) <- names(past_berpar) <- rownames(past_results)
  expect_equal(round(margalef(data_agr, tax_lev = "Taxa"), 3),  past_mar)
  expect_equal(round(brill(data_agr, tax_lev = "Taxa"), 3),  past_brill)
  expect_equal(round(menhinick(data_agr, tax_lev = "Taxa"), 3),  round(past_men, 3) )
  expect_equal(round(berpar(data_agr, tax_lev = "Taxa"), 3),  round(past_berpar, 3))
  expect_equal(round(invberpar(data_agr, tax_lev = "Taxa"), 2),  round(1/past_berpar, 2))
  expect_equal(round(invberpar(data_agr, tax_lev = "Taxa"), 2),  round(1/past_berpar, 2))
})


