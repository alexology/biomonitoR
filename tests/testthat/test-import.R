test_that("unassigned", {
  data( macro_ex )
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

