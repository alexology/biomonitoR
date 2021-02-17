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


# we should test genus results from past too

test_that("past_abu", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  # results compared with those of the past software version 4.04 accessed on 2021-01-11
  load(system.file("testdata", "past_results_abu.rda", package="biomonitoR"))
  past_results <- as.data.frame (past_results_abu[,-1], row.names = past_results_abu[,1])
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
})

test_that("past_bin", {
  data(macro_ex)
  data_bio <- suppressWarnings(as_biomonitor(macro_ex, FUN = bin))
  data_agr <- aggregate_taxa(data_bio)
  # results compared with those of the past software version 4.04 accessed on 2021-01-11
  load(system.file("testdata", "past_results_bin.rda", package="biomonitoR"))
  past_results <- as.data.frame (past_results_bin[,-1], row.names = past_results_bin[,1])
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
})

test_that("past_genus", {
  data(macro_ex)
  data_bio <- suppressWarnings(as_biomonitor(macro_ex, FUN = sum))
  data_agr <- aggregate_taxa(data_bio)
  # results compared with those of the past software version 4.04 accessed on 2021-01-11
  load(system.file("testdata", "past_results_genus.rda", package="biomonitoR"))
  past_results <- as.data.frame (past_results_genus[,-1], row.names = past_results_genus[,1])
  past_mar <- past_results[, "Margalef"]
  past_brill <- past_results[, "Brillouin"]
  past_men <- past_results[, "Menhinick"]
  past_berpar <- past_results[, "Berger-Parker"]
  names(past_mar) <- names(past_brill) <- names(past_men) <- names(past_berpar) <- rownames(past_results)
  expect_equal(round(margalef(data_agr, tax_lev = "Genus"), 3),  past_mar)
  expect_equal(round(brill(data_agr, tax_lev = "Genus"), 3),  past_brill)
  expect_equal(round(menhinick(data_agr, tax_lev = "Genus"), 2),  round(past_men, 2) )
  expect_equal(round(berpar(data_agr, tax_lev = "Genus"), 3),  round(past_berpar, 3))
  expect_equal(round(invberpar(data_agr, tax_lev = "Genus"), 3),  round(1/past_berpar, 3))
})



test_that("abdiv", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  apply(macro_ex[, -1], 2, function(x) abdiv::mcintosh_d(x))
  expect_equal(mcintosh(data_agr, tax_lev = "Taxa"), apply(macro_ex[, -1], 2, function(x) abdiv::mcintosh_d(x)))
  data_veg_gen <- convert_to_vegan(data_agr, tax_lev = "Genus")
  expect_equal(mcintosh(data_agr, tax_lev = "Genus"), apply(data_veg_gen, 1, function(x) abdiv::mcintosh_d(x)))
  data_bio_bin <- suppressWarnings(as_biomonitor(macro_ex, FUN = bin))
  data_agr_bin <- aggregate_taxa(data_bio_bin)
  data_veg_gen_bin <- convert_to_vegan(data_agr_bin, tax_lev = "Genus")
  expect_equal(mcintosh(data_agr_bin, tax_lev = "Genus"), apply(data_veg_gen_bin, 1, function(x) abdiv::mcintosh_d(x)))
})


test_that("allindices", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)

  data.sha <- shannon(data_agr)
  data.berpar <- berpar(data_agr)
  data.brill <- brill(data_agr)
  data.invb <- invberpar(data_agr)
  data.invs <- invsimpson(data_agr)
  data.mar <- margalef(data_agr)
  data.mch <- mcintosh(data_agr)
  data.men <- menhinick(data_agr)
  data.pie <- pielou(data_agr)
  data.sim <- simpson(data_agr)
  data.esi <- esimpson(data_agr)
  data.fis <- fisher(data_agr)

  res <- data.frame(data.sha, data.berpar, data.brill, data.invb, data.invs, data.mar, data.mch,
                    data.men, data.pie, data.sim, data.esi, data.fis)

  colnames(res) <- c(
    "shannon", "berpar", "brill", "invberpar", "invsimpson",
    "margalef", "mcintosh", "menhinick", "pielou",
    "simpson", "esimpson", "fisher"
  )

  expect_equal(res, allindices(data_agr))

})

