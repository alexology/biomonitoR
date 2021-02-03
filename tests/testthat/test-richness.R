test_that("basic_richness", {
  data( macro_ex )
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_genus <- c(13,20)
  names( data_genus ) <- names( macro_ex )[ -1 ]
  expect_equal( suppressWarnings( richness( data_agr , tax_lev = "Genus" ) ) ,  data_genus )
  expect_equal( richness(data_agr , tax_lev = "Taxa") ,  apply( macro_ex[ , -1 ] , 2 , function( x ) sum( x > 0 ) ) )
  expect_equal( abundance(data_agr , tax_lev = "Taxa") ,  apply( macro_ex[ , -1 ] , 2 , sum ) )
})


test_that( "get_taxa_richness", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_genus <- c(13,20)
  names(data_genus) <- names(macro_ex)[ -1 ]
  expect_equal(get_taxa_abundance( data_agr , taxa = "Ancylus" ),
                t( macro_ex[ macro_ex$Taxa %in% "Ancylus" , -1 ] )[ , 1 ]  )

  ancy_abu <- t( macro_ex[ macro_ex$Taxa %in% "Ancylus" , -1 ] )[ , 1 ]
  expect_equal(get_taxa_abundance( data_agr , taxa = "Ancylus", rel = TRUE), ancy_abu /apply(macro_ex[,-1],2,sum))
  expect_error(get_taxa_abundance(data_agr), "Please provide a taxon name")
  expect_error(get_taxa_abundance(data_agr, taxa = "Ergo"), "None of the taxa provided were found in the data_agr database")


  macro_ex_bin <- to_bin(macro_ex)
  data_bio_bin <- as_biomonitor(macro_ex_bin, FUN = bin)
  data_agr_bin <- aggregate_taxa(data_bio_bin)
  expect_error(get_taxa_abundance( data_agr_bin , taxa = "Ancylus" ), "This function cannot be applied to presence-absence data.")

})
