#' psi
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates the *Proportion of Sediment-sensitive Invertebrates index* (PSI) according to the most recent version used in UK.
#' @param x result of the function aggregatoR.
#' @param method the only choice is `extence`. Users can provide their own data.frame (see examples) with a column called *Taxon* and the column of scores called *Score*.
#' @param abucl Log abundance categories. Treshold are set to 1, 9, 99 and 999 as in the original paper of Extence et al. (2013).
#' @param agg this option allows the composite family approach. It can be `FALSE`, `TRUE` or a `data.frame`.
#' If `FALSE` no aggregation will be performed, while if `TRUE` aggregation will be performed according to the rules described in Details.
#' A `data.frame` containing the aggregation rules can be provided by the user.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' `agg` cannot be `TRUE` when a `data.frame` is provided as method.
#' @param exceptions taxa that need to be exluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' `agg` cannot be `TRUE` when a data.frame is provided as `method`.
#' @param fssr_scores Optional, scores (fssr) for different abundance categories of taxa associated with
#'  *Fine Sediment Sensitivity Ratings*. To be used when a custom `method` is provided.
#' @param traceB if set to `TRUE` a list as specified below will be returned.
#'
#' @keywords psi
#' @details Despite Extence et al 2013 did not suggest nay aggregation rule, the `psi` implementation of `biomonitoR`
#' allows for a default aggregation as specified below. Custom aggregation rules can be provided as a `data.frame`.
#' \enumerate{
#'  \item Tipulidae (inc. Limoniidae, Pediciidae & Cylindrotomidae)
#'  \item Siphlonuridae (inc. Ameletidae)
#'  \item Hydrophilidae (inc. Georissidae, Helophoridae & Hydrochidae)
#'  }
#'
#' The `psi` function automatically check for parent-child pairs in the scoring system, see the return section for a definition.
#' All the information used for `psi` calculation can be retrieved with the function \code{\link{showscores}}.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of the `psi` index.
#'  \item `taxa_df` The data.frame used for the calculation containing the abundance of taxa receiving a score.
#'  \item `abu_df` The data.frame containing fssr scores and abundance classes for each site.
#'  \item `psi_df` The data.frame used for the calculation containing scores for each site.
#'  \item `composite_taxa` Taxa aggregated following the aggregation rules when agg is not `NULL`.
#'  \item `exceptions` A data.frame containing the con.
#'  \item `parent_child_pairs` For instance in Spanish `aspt` both Ferrissia and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#'
#' @references Extence CA, Chadd RP, England J, Dunbar MJ, Wood PJ, Taylor ED. 2013. The assessment of fine sediment accumulation in rivers using macro-invertebrate community response. River Research and Applications 29, 17-55.
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated PSI scores and algorithms.
#' @importFrom dplyr '%>%' select inner_join group_by summarise rename filter
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble deframe
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' psi( data.agR )
#'
#' # change abundance classes
#' psi( data.agR , abucl = c( 1 , 9 , 99 , 999 , 9999 ) )
#'
#' # provide your own score system. Scores and aggregation rules are for example purpose only.
#'
#' psi_fssr <- data.frame( Taxon = c( "Ephemerellidae" , "Leuctridae" , "Chironomidae" ) ,
#'  FSSR_Score = c( 1 , 2 , 3 ) )
#'
#' psi_acc <- data.frame( Taxon = "Ephemerellidae" , Correct_Taxon = "Chironomidae" )
#'
#' fssr_scores <- data.frame( FSSR = rep( 1:3 , each = 3 ) , ABUCLASS = rep( 1:3 , 3 ) ,
#' SCORE = c( 9 , 10 , 11 , 8 , 9 , 10 , 7 , 7 , 7 ) )
#'
#' # without aggregation rules
#' psi( data.agR , method = psi_fssr , fssr_scores = fssr_scores , traceB = TRUE )
#'
#' # with aggregation
#'
#' psi( data.agR , method = psi_fssr , agg = psi_acc , fssr_scores = fssr_scores , traceB = TRUE )
#'

psi <- function( x , method = "extence" , abucl = c( 1 , 9 , 99 , 999 ) , agg = FALSE , fssr_scores = NULL , exceptions = NULL , traceB = FALSE ){


  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # useful for transforming data to 0-1 later
  if( inherits( x , "bin" ) ){
    BIN <- TRUE
  } else { BIN <- FALSE }


  if( ! any( identical( method , "extence" ) | is.data.frame( method ) ) ) stop( "method is not extence or a custom data.frame" )
  if( ! any( isFALSE( agg ) | isTRUE( agg ) | is.data.frame( agg ) ) ) stop( "agg is not one of TRUE, FALSE or a custom data.frame" )
  if( ! any( is.null( fssr_scores ) | is.data.frame( fssr_scores ) ) ) stop( "fssr_scores is not one of NULL or a custom data.frame" )

  # Store tree for searching for inconsistencies
  Tree <- x[[ "Tree" ]][ , 1:10 ]

  numb <- c( which( names( x ) == "Tree" ) , which( names( x ) == "Taxa" ) ) # position of the Tree and Taxa data.frame in the biomonitoR object that need to be removed

  # remove Tree and Taxa data.frame
  x <- x[ -numb ]
  st.names <- names( x[[ 1 ]][ -1 ] ) # names of the sampled sites

  # create dummy variables to avoid R CMD check NOTES


  z <- w <- Taxon <- Score <- Sample <- FSSR <- FSSR_SCORE <- FSSR_Score <- SCORE <-   NULL


  # the following if statement is to allow the users to provide their own psi scores and aggregation rules.
  # y represents the method to be used

  if( is.data.frame( method ) ){
    if( ! ( isFALSE( agg ) | is.data.frame( agg ) ) ){
      stop( "When method is a data.frame agg needs to be FALSE or a data.frame containing the aggregation rules")
    }
    y <- method
    if( ! is.data.frame( fssr_scores ) ) stop( "fssr_scores needed when method is a data.frame provided by the user" )
    w <- fssr_scores
    if( is.data.frame( agg ) ){
      z <- agg
    }

  } else {

    if( ! ( isTRUE( agg ) | isFALSE( agg ) ) ) stop( "When using the deafult method agg can only be TRUE or FALSE")
    if( ! is.null( fssr_scores ) ) stop( "When using deafult methods fssr_scores can only be NULL")

    # assign the default scores and aggregation rules as needed by the user
    if( identical( method , "extence" ) ){
      y <- psi_scores_fam_extence
      w <- psi_fssr
      if( isTRUE( agg ) ){
        z <- epsi_acc_fam_uk
      }
    }
  }


  # the calculation of the index in biomonitoR consists in rbind all the taxonomic levels
  # in the biomonitoR object that has been previously deprived of Taxa and Tree elements and then merge
  # it with the scores data.frame.
  # The first step is to change the column name of the first column of each data.frame to
  # an unique name

  for( i in 1:length( x ) ){
    colnames( x[[ i ]] )[ 1 ] <- "Taxon"
  }

  # rbind the data.frames representing a taxonomic level each
  # aggregate is not necessary here
  DF <- do.call( "rbind" , x )
  rownames( DF ) <- NULL
  DF <- aggregate(. ~ Taxon, DF , sum)

  if( ! is.null( exceptions ) ){
    DF <- manage_exceptions( DF = DF , Tree = Tree , y = y , Taxon = exceptions )
    if( ! is.data.frame( DF ) ){
      exce <- DF[[ 2 ]]
      DF <- DF[[ 1 ]]
    }
  }

  # merge the new data.frame with the score data.frame and change
  # the names of the taxa according to the aggregation rules if needed
  DF <- merge( DF, y[ , "Taxon" , drop = FALSE ] )
  if( ! is.null( z ) ){
    taxa.to.change <- as.character( DF$Taxon[ DF$Taxon %in% z$Taxon  ] )
    DF <- checkBmwpFam( DF = DF , famNames = z , stNames = st.names )
  } else {
    DF <- DF
  }

  if( traceB ){
    df2 <- DF
  }

  DF <- manage_inconsistencies( DF = DF , Tree = Tree )
  if( ! is.data.frame( DF ) ){
    incon <- DF[[ 2 ]]
    DF <- DF[[ 1 ]]
  }

  # transform the data.frame from abundance to presence-absence if needed
  if( BIN ){
    DF <- to_bin( DF )
  }

  DF <- merge( y , DF )

  class.fun <- function( x ) cut( x , breaks = c( abucl , 10^18 ) , labels = 1:length( abucl ) , include.lowest = TRUE , right = TRUE )
  abu.class <- apply( apply( DF[ , -c( 1:2 ) , drop = FALSE ] , 2 , class.fun ) , 2 , as.numeric )
  abu.class[ is.na( abu.class ) ] <- 0
  tot.mer <- data.frame( DF[ , 1:2 ] , abu.class , check.names = FALSE )


  res.tot <- tot.mer %>% pivot_longer( -c( Taxon , FSSR_Score ) , names_to = "Sample", values_to = "ABUCLASS" ) %>%
    dplyr::rename( "FSSR" = "FSSR_Score" ) %>% inner_join( w , by = c( "FSSR" , "ABUCLASS" ) ) %>% select( c( "Sample" , "SCORE" ) ) %>%
    group_by( Sample ) %>% summarise( PSI = sum( SCORE ) ) %>% deframe()

  res.sen <- tot.mer %>% pivot_longer( -c( Taxon , FSSR_Score ) , names_to = "Sample", values_to = "ABUCLASS" ) %>%
    dplyr::filter( FSSR_Score == 1 | FSSR_Score == 2 ) %>%
    dplyr::rename( "FSSR" = "FSSR_Score" ) %>% inner_join( w , by = c( "FSSR" , "ABUCLASS" ) ) %>% select( c( "Sample" , "SCORE" ) ) %>%
    group_by( Sample ) %>% summarise( PSI = sum( SCORE ) ) %>% deframe()

  if( length( res.sen ) == 0 ) ( res.sen <- rep( 0 , length( st.names ) ) )

  # assure to have results if all indicator taxa are missing from a sample
  if( length( res.tot ) != length( st.names ) | length( res.sen ) != length( st.names ) ){
    temp <- merge( data.frame( Sites = st.names ) , data.frame( Sites = names( res.tot ) , res_tot = res.tot ) , all = TRUE )
    temp <- merge( temp ,  data.frame( Sites = names( res.sen ) , res_sen = res.sen ) , all = TRUE )
    res <- temp$res_sen / temp$res_tot * 100
    names( res ) <- temp$Sites
  } else {
    res <- res.sen / res.tot * 100
    names( res ) <- st.names
  }


  res <- res[ match( st.names , names( res ) ) ]

  if( traceB == FALSE ){
    res
  } else {
    if( ! exists( "taxa.to.change" , inherits = FALSE ) ){
      df3 <- "none"
    } else { df3 <- taxa.to.change }
    if( exists( "exce" , inherits = FALSE  ) ){
      df4 <- exce
    } else { df4 <- "none" }
    if( exists( "incon" , inherits = FALSE  ) ){
      df5 <- incon
    } else { df5 <- "none" }

    temp <- tot.mer %>% pivot_longer( -c( Taxon , FSSR_Score ) , names_to = "Sample", values_to = "ABUCLASS" ) %>%
      dplyr::rename( "FSSR" = "FSSR_Score" ) %>% inner_join( w , by = c( "FSSR" , "ABUCLASS" ) ) %>% select( c( "Taxon" , "Sample" , "SCORE" ) ) %>%
      pivot_wider( names_from = Sample , values_from = SCORE , values_fill = list( SCORE = 0 ) ) %>% as.data.frame( )
    res <- list( results = res , taxa_df = df2 , abu_cl = tot.mer , psi_df = temp , composite_taxa = df3 ,  exceptions = df4 ,  parent_child_pairs = df5 )
    res
  }
}
