#' traitScaling
#'
#' A function for scaling traits across taxonomic levels.
#'
#' This function allows to obtain missing traits for a target taxon by taking traits from lower or to upper taxomic levels.
#' For instance, consider the case where the genus Acroloxus is present in the user dataset and the species Acroloxus lacustris
#' in the traits database. A simple merge would exclude Acroloxus from the tha analysis since Acroloxus and A. lacustris
#' would not match. The function traitScaling allows to assign Acroloxus lacustris traits to Acroloxus.
#' This function works also in the opposite direction. Consider the case where there are no traits for the target taxon and
#' the target taxon has not been identified at species level. The function traitScaling will assign the traits of the nearest taxonomic level
#' to the target taxa (e.g. Tanypodinae traits assign to Ablabesmyia monilis). Consider also these examples to understand the behaviour of this
#' function. For instance Anabolia lombarda is present in the user taxomic dataset while only Anabolia nervosa and Anabolia are present
#' in the trait database. In this case traitScaling will assing to A. lombarda only the traits of Anabolia.
#' Moreover, let assume that Coelostoma is present in the user dataset while only Berosus and Crenitis punctatostriata are present in the traits database.
#' Here traitScaling will assign to Coelostoma the scores of Berosus and C. punctatostriata since they belong to the same family and there are no information at family level.\cr
#' In case of multiple matches for the target taxa, traitScaling will not average the scores but it will provide all the matches. It is
#' then straightforward to average the results by using the \link[stats]{aggregate}. \cr
#' The function traitScaling will measure the taxonomic distance between the target taxa and the taxa used to assign the trait score. This distance
#' can be negative (e.g. Species to Genus) and positive (e.g. Genus to Species). The distance is measured assigning values as follows:
#' 1 (Species to Genus) , 2 (Species to family), -1 (Family to Genus), etc. traitScaling considers only the
#' taxonomic levels from Subspecies to Family (Subspecies, Species, Genus, Tribus, Subfamily, Family).
#'
#'
#' @param x results of function aggregatoR
#' @param traitDB a trait data base with a column `Taxa` and the other columns
#'   containing the traits. If a trait has several modalities they should be
#'   named as follow: TRAIT_MODALITY.
#'
#'   By default, the data base used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#' @param dfref reference database as used in the function aggregatoR.
#' @param filter_by_distance filter the results according to the taxonomic distance. Possible values are "pos" , "neg" or a positive integer. See details.
#' @note USE WITH CAUTION, STILL IN DEVELOPMENT.
#'
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.ts <- traitScaling( data.agR )
#'
#' # averaging
#'
#' data.ts.av <- aggregate( . ~ Taxa , data = data.ts[ , -c(2:5)] , FUN = mean  )
#'
#' # traits random sampling
#' data.ts.st <- sampleTraits( data.ts )
#' @seealso [aggregatoR]
#' @export
#' @export sampleTraits



traitScaling <-  function( x , traitDB = NULL , dfref = NULL , filter_by_distance = NULL ){

  if( is.null( traitDB ) ){
    # check if x is of class biomonitoR and mi
    classCheck( x , group = "mi")

    traitDB <- traitsTachet

    # create dummy variables to avoid R CMD check NOTES
    traitsTachet <- Taxa <- modality <- affinity <- Phylum <- Subspecies <-
      Abundance <- Sample <- Weight <- Affinity <- totWeight <-
      weightedAffinity <- Category <- . <- NULL

    # prepare the taxa trait database
    trait_db <- traitDB                               %>%
      (function(df) {
        mutate(df,
               Taxa = gsub(pattern     = "sp[.]|Ad[.]|Lv[.]|Gen[.]|lv[.]|ad[.]|gen[.]",
                           replacement = "",
                           x           = Taxa))
      })                                              %>%
      gather(key = modality, value = affinity, -Taxa) %>%
      group_by(Taxa, modality)                        %>%
      summarise(affinity = mean(affinity))            %>%
      spread(key = modality, value = affinity)        %>%
      ungroup()
    trait_db$Taxa <- trimws(trait_db$Taxa)
    trait_db <- as.data.frame( trait_db )
    trait_db <- trait_db[ , match( colnames(traitDB) , colnames( trait_db) ) ]
  } else{
    trait_db <- traitDB
    trait_db$Taxa <- trimws(trait_db$Taxa)
    classCheck(x)
  }

  if( is.null( dfref ) ){
    dfref <- mi_ref
  } else { dfref <- dfref }



  # merge the trait database with the reference database in order to scale
  # traits across taxonomic levels
  # [ , - 1] deletes the Taxa columns
  ref <- merge( dfref , trait_db , by = "Taxa" , sort = FALSE )[ , -1 ]

  # create a data.frame with the same column as trait_db but with 0 rows
  # it will be important later to iterate rbind to this object
  trait.interm <- data.frame( trait_db[-c(1:nrow(trait_db)) , ] )

  # ref.na is a 0 length vector to store the taxa names of the selected rows
  # this allow to does not loose the information about the name of the taxa at the orginial
  # taxonomic level

  ref.na <- c()

  # 10 to  5  because it is intended to work from subspecises to family

  # cycle to scale the traits among taxonomic levels
  for( i in 10:5){
    temp <- ref[ , -which( c(1:10) != i ) ]
    names( temp )[ 1 ] <- "Taxa"
    temp <- temp[ temp[ , 1 ] != "" , ]
    ref.name <- ref[ rownames( ref ) %in% rownames( temp ) , 5:10 ]
    ref.name <- apply( ref.name , 1 , function( x )(  rev( x )[ rev( x ) != "" ][ 1 ]  ) )
    ref.na <- c( ref.na , ref.name )
    trait.interm <- rbind( trait.interm , temp)
  }

  trait_db <- trait.interm

  DF <- x[["Tree"]]

  # DFtaxa stores the taxa present in the user database

  DFtaxa <- as.character( DF[  11 ] )
  result.list <- apply( DF , 1 , function( x ) traitS( x = x, y = DF , z = trait_db , w = ref.na ) )
  result.data.frame <- do.call( rbind , result.list )
  result.data.frame.single <- result.data.frame

  unique.taxa <- unique( result.data.frame.single$Taxa_db )

  # deleting row for the follwing reason. It happens that traits are present for a species
  # and also for the genus of this species. If the user sample
  # contains a species other than that reported in the trait list we want only to keep the
  # trait at genus level

  for( i in 1:length( unique.taxa ) ){
    res <- result.data.frame.single[ result.data.frame.single$Taxa_db %in% unique.taxa[ i ] , 1:3 ]
    res.sum <- ( res[ , 1 ] == res[ , 2 ] ) + ( res[ , 1 ] == res[ , 3 ] ) + ( res[ , 2 ] == res[ , 3 ] )
    if( any( res.sum == 0 ) ){
      if( sum( res.sum )!= 0){
        to.del <- which( result.data.frame.single$Taxa_db %in% unique.taxa[ i ] )[ res.sum == 0 ]
        result.data.frame.single <- result.data.frame.single[ -to.del , ]
      }
    }
  }

  ref_long <- data.frame( Taxonomic_level = character( ) , Taxa = character( ) , stringsAsFactors = FALSE )

  for( i in 10:1 ){
    temp <- as.character( dfref[ , i ] )
    temp <- temp[ temp != ""]
    temp.rep <- rep( names(ref[ , i , drop = FALSE] ) , length( temp ) )
    temp.df <- data.frame( Taxonomic_level = temp.rep , Taxa = temp )
    ref_long <- rbind( ref_long , temp.df )
  }

  ref_long <- ref_long[ !duplicated( ref_long ) , ]

  taxa_db.taxlev <- ref_long[ match( ref_long$Taxa , result.data.frame.single$Taxa_db ) ,  ]

  result.data.frame.single$Taxa_db <- as.character( result.data.frame.single$Taxa_db )
  result.data.frame.single$Traits_real <- as.character( result.data.frame.single$Traits_real )
  ref_long$Taxa <- as.character( ref_long$Taxa )

  taxa_db.taxlev <- inner_join( result.data.frame.single[ , 1 , drop = FALSE ] , ref_long , by = c( "Taxa_db" = "Taxa"))
  traits.taxlev <- inner_join( result.data.frame.single[ , 2 , drop = FALSE ] , ref_long , by = c( "Traits_real" = "Taxa"))
  names( taxa_db.taxlev ) <- c( "Taxa_taxlev" , "Taxa_db" )
  names( traits.taxlev ) <- c("Traits_taxlev" , "Traits_real" )


  taxLev.info <- data.frame( taxa_db.taxlev , traits.taxlev , stringsAsFactors = FALSE )

  dist.taxlev <- data.frame( taxLev = names( dfref )[ -11 ] , distance = c(10:1) )
  dist.taxlev$taxLev <- as.character( dist.taxlev$taxLev )
  taxLev.info$Taxa_db <- as.character( taxLev.info$Taxa_db )
  taxLev.info$Traits_real <- as.character( taxLev.info$Traits_real )

  a <- inner_join( taxLev.info[ , c( 2 , 4 ) ] , dist.taxlev , by = c( "Taxa_db" = "taxLev"))[ , 3 ]
  b <- inner_join( taxLev.info[ , c( 2 , 4 ) ] , dist.taxlev , by = c( "Traits_real" = "taxLev"))[ , 3 ]

  taxLev.info$Taxonomic_distance <- a-b

  names( taxLev.info )[ 1 ] <- "Taxa"
  final.traits <- data.frame( taxLev.info , result.data.frame.single[ , -c(1:3) ])
  rownames( final.traits ) <- NULL
  if( is.null( filter_by_distance ) ){
    final.traits
  } else {
    if( is.character( filter_by_distance) ){
      if( filter_by_distance == "pos" ){ final.traits[ final.traits$Taxonomic_distance >= 0, ] } else { if( filter_by_distance == "neg" )
      { final.traits[ final.traits$Taxonomic_distance <= 0 , ] } else {
        stop("pos, neg or an integer are needed when filter_by_distance is not NULL") }
      }
    } else {
      final.traits[ final.traits$Taxonomic_distance <= abs( filter_by_distance ) , ]
    }
  }
}

