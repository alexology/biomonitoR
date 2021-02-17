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


test_that( "get_taxa_abundance", {
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



test_that( "get_taxa_richness", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)

  data_tree <- to_bin(data_agr[["Tree"]])
  data_tree_ephe <- subset(data_tree, Order == "Ephemeroptera")
  data_tree_genus <- aggregate(. ~ Genus, data_tree_ephe[ , c("Genus", "Sample_1", "Sample_2")], FUN = bin)
  data_tree_family <- aggregate(. ~ Family, data_tree_ephe[ , c("Family", "Sample_1", "Sample_2")], FUN = bin)
  get_taxa_richness(data_agr, taxa = "Diptera", tax_lev = "Subfamily")


  expect_equal(get_taxa_richness(data_agr, taxa = "Ephemeroptera", tax_lev = "Taxa"), apply(data_tree[ data_tree[, "Order"] %in% "Ephemeroptera" , ][, -c(1:11)], 2, sum))
  expect_equal(get_taxa_richness(data_agr, taxa = "Ephemeroptera", tax_lev = "Genus"), apply(data_tree_genus[, -1], 2, sum))
  expect_equal(get_taxa_richness(data_agr, taxa = "Ephemeroptera", tax_lev = "Family"), apply(data_tree_family[, -1], 2, sum))
  expect_equal(get_taxa_richness(data_agr, taxa = c("Baetidae", "Caenidae", "Ephemerellidae", "Heptageniidae", "Leptophlebiidae"), tax_lev = "Family"), apply(data_tree_family[, -1], 2, sum))
  expect_error(get_taxa_richness(data_agr, tax_lev = "Family"), "Please provide a taxon name and/or a taxonomic level")
  expect_error(get_taxa_richness(data_agr), "Please provide a taxon name and/or a taxonomic level")
  expect_error(get_taxa_richness(data_agr, taxa = c("Ephemeroptera", "Trichoptera", "Diptera"), tax_lev = c("Family", "Genus")), "tax_lev must be of the same length of taxa")
  expect_error(suppressWarnings(get_taxa_richness(data_agr, taxa = "Baetis", tax_lev = "Family")), "Taxonomic level of taxa cannot be lower than taxonomic level of tax_lev")
  expect_error(get_taxa_richness(data_agr, taxa = "Epemeroptera", tax_lev = "Family"), "Please provide a valid taxon name. Names provided can also be absent in your database.")

})



test_that( "all_rich", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)

  data_allrich <- suppressWarnings(allrich(data_agr))
  data_phy <- suppressWarnings(richness(data_agr, tax_lev = "Phylum"))
  data_cla <- suppressWarnings(richness(data_agr, tax_lev = "Class"))
  data_scl <- suppressWarnings(richness(data_agr, tax_lev = "Subclass"))
  data_ord <- suppressWarnings(richness(data_agr, tax_lev = "Order"))
  data_fam <- suppressWarnings(richness(data_agr, tax_lev = "Family"))
  data_sfa <- suppressWarnings(richness(data_agr, tax_lev = "Subfamily"))
  data_tri <- suppressWarnings(richness(data_agr, tax_lev = "Tribus"))
  data_gen <- suppressWarnings(richness(data_agr, tax_lev = "Genus"))
  data_spe <- suppressWarnings(richness(data_agr, tax_lev = "Species"))
  data_ssp <- suppressWarnings(richness(data_agr, tax_lev = "Subspecies"))
  data_tax <- suppressWarnings(richness(data_agr, tax_lev = "Taxa"))

  res <- data.frame(Phylum = data_phy, Class = data_cla, Subclass = data_scl, Order = data_ord,
                    Family = data_fam, Subfamily = data_sfa, Tribus = data_tri, Genus = data_gen,
                    Species = data_spe, Subspecies = data_ssp, Taxa = data_tax)

  expect_equal(data_allrich, res)
})
