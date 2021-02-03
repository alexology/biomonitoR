test_that("unassigned", {
  data(macro_ex)
  load(system.file("testdata", "ref_def.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  macro_ex_bin <- to_bin(macro_ex)
  data_bio_bin <- suppressMessages(as_biomonitor(macro_ex_bin, FUN = bin))
  data_agr_bin <- aggregate_taxa(data_bio_bin)
  data_bio_custom <- as_biomonitor(macro_ex, dfref = ref_def)
  data_agr_custom <- aggregate_taxa(data_bio_custom)
  data_agr_genus <- data_agr[["Genus"]][-1, ]
  names(data_agr_genus)[1] <- "Taxa"
  rownames(data_agr_genus) <- NULL
  expect_equal(data_agr[["Genus"]][ 1 , 1 ] ,  "unassigned" )
  expect_equal(length( data_agr ) ,  12 )
  expect_equal(class( data_bio ) ,  "asb" )
  expect_equal(class( data_agr ) ,  "biomonitoR" )
  expect_equal(names( data_agr )[1:11] , c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa") )
  expect_equal(class( data_agr_bin ) ,  c("biomonitoR","bin"))
  expect_equal(any(macro_ex_bin[,-1]>1),  FALSE)
  expect_equal(any(as.data.frame(data_bio_bin)[,-c(1:11)]>1),  FALSE)
  expect_warning(as_biomonitor(macro_ex_bin, FUN = sum), "Presence-absence data detected but FUN is not set to bin. Is it this what you want?")
  expect_equal(class( data_agr_custom ) ,  c("biomonitoR","custom"))
  load(system.file("testdata", "macro_ex_genus.rda", package="biomonitoR"))
  expect_equal(data_agr_genus,  macro_ex_genus)
})


test_that("to_change", {
  data( macro_ex )
  load(system.file("testdata", "to_change_asb.rda", package="biomonitoR"))
  load(system.file("testdata", "macro_ex_mod_asb.rda", package="biomonitoR"))
  data_bio <- as_biomonitor(macro_ex, to_change = to_change_asb)
  data_agr <- aggregate_taxa(data_bio)
  data_agr <- data_agr[["Taxa"]]
  names(data_agr)[1] <- "Taxa"
  expect_equal(data_agr,  macro_ex_mod_asb)
})


test_that("checks", {
  data(macro_ex)
  macro_ex_fail1 <- macro_ex_fail2 <- macro_ex_NA <-  macro_ex
  names(macro_ex_fail1)[1] <- "Fail"
  macro_ex_fail2[, 2] <- as.factor(macro_ex_fail2[, 2])
  to_change_fail1 <- data.frame(Taxon = c("Baetis", "Baetis"), Correct_Taxon = c("Fail1", "Fail2"))
  to_change_fail2 <- data.frame(Taxon = character(), Correct_Taxon = character())
  macro_ex_NA[macro_ex_NA == 0] <- NA
  expect_error(as_biomonitor(macro_ex_fail1), "A column called Taxa is needed")
  expect_error(as_biomonitor(macro_ex_fail2), "Non-numeric columns other than Taxa are not allowed")
  expect_error(as_biomonitor(macro_ex, to_change = to_change_fail1), "the same name cannot be present twice in the Taxon column of the data.frame to_change")
  expect_error(as_biomonitor(macro_ex, to_change = to_change_fail2), "to_change must have at least one entry")
  expect_error(as_biomonitor(macro_ex, to_change = c(1:3)), "to_change needs to be NULL, default or data.frame as specified in the help")
  expect_message(as_biomonitor(macro_ex_NA), "NA detected, transformed to 0")
})


test_that("default_trace_back", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex, traceB = TRUE)


  expect_equal(length(data_bio),  1)
  expect_equal(names(data_bio),  "taxa_db")

  to_change_fake <- data.frame(Taxon = "Acentrella", Correct_Taxon = "Acentrella")
  data_bio_f <- as_biomonitor(macro_ex, to_change = to_change_fake, traceB = TRUE)
  expect_equal(length(data_bio_f),  2)
  expect_equal(names(data_bio_f),  c("taxa_db", "corrected_names"))


})



test_that("trace_back", {
  load(system.file("testdata", "metrics_uk_data.rda", package="biomonitoR"))
  load(system.file("testdata", "taxa_uk_data.rda", package="biomonitoR"))
  load(system.file("testdata", "to_change.rda", package="biomonitoR"))

  expect_message(as_biomonitor(taxa_uk_data), "Some taxa were excluded, check with traceB = TRUE")


  taxa_uk_data_asb <- suppressMessages(as_biomonitor(taxa_uk_data, to_change = NULL, traceB = TRUE))
  expect_equal(length(taxa_uk_data_asb),  2)
  expect_equal(names(taxa_uk_data_asb),  c("taxa_db", "suggested_taxa_names"))

  taxa_uk_data_asb_tc1 <- suppressMessages(as_biomonitor(taxa_uk_data, to_change = "default", traceB = TRUE))
  expect_equal(length(taxa_uk_data_asb_tc1),  3)
  expect_equal(names(taxa_uk_data_asb_tc1),  c("taxa_db", "corrected_names", "suggested_taxa_names"))

  taxa_uk_data_asb_tc2 <- suppressMessages(as_biomonitor(taxa_uk_data, to_change = to_change, traceB = TRUE))
  expect_equal(length(taxa_uk_data_asb_tc2),  3)
  expect_equal(names(taxa_uk_data_asb_tc2),  c("taxa_db", "corrected_names", "suggested_taxa_names"))


})



test_that("macrophytes", {
  data(oglio)
  oglio_asb <- suppressWarnings(as_biomonitor(oglio, group = "mf"))
  oglio_df <- as.data.frame(oglio_asb)
  expect_equal(nrow(oglio),  nrow(oglio_df))

})


test_that("fish", {
  fish_data <- data.frame(Taxa = c("Cottus gobio", "Silurus glanis"), site_1 = c(1,0), site_2 = c(1,1))
  fish_asb <- suppressWarnings(as_biomonitor(fish_data, group = "fi", FUN = bin))
  fish_df <- as.data.frame(fish_asb)
  expect_equal(nrow(fish_data),  nrow(fish_df))

})


test_that("class_aggregate_taxa", {
  expect_error(aggregate_taxa(c(1:3)), "x is not an object created with as_biomonitor")
})
