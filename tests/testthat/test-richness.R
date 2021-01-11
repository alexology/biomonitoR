test_that("basic_richness", {
  data( macro_ex )
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_genus <- c(13,20)
  names( data_genus ) <- names( macro_ex )[ -1 ]
  expect_equal( suppressWarnings( richness( data_agr , tax_lev = "Genus" ) ) ,  data_genus )
  expect_equal( richness( data_agr , tax_lev = "Taxa" ) ,  apply( macro_ex[ , -1 ] , 2 , function( x ) sum( x > 0 ) ) )
  expect_equal( abundance( data_agr , tax_lev = "Taxa" ) ,  apply( macro_ex[ , -1 ] , 2 , sum ) )
})


test_that( "get_taxa_richness", {
  data( macro_ex )
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_genus <- c(13,20)
  names( data_genus ) <- names( macro_ex )[ -1 ]
  expect_equal( get_taxa_abundance( data_agr , taxa = "Ancylus" ) ,
                t( macro_ex[ macro_ex$Taxa %in% "Ancylus" , -1 ] )[ , 1 ]  )
})
