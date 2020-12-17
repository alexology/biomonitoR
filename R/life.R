#' life
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates LIFE index according to most recent version used in UK.
#'
#' @param x result of the function aggregatoR.
#' @param method possible choices are "extence" and "life_2017". A custom `data.frame` containing the scoring system can be provided by the user.
#' @param abucl Log abundance categories. Treshold are set to 1, 9, 99, 999 and 9999 as in the original paper of Extence et al. (1999).
#' @param agg this option allows the composite family approach. It can be `FALSE`, `TRUE` or a `data.frame`.
#' If `FALSE` no aggregation will be performed, while if `TRUE` aggregation will be performed according to the rules described in Details.
#' A `data.frame` containing the aggregation rules can be provided by the user.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' `agg` cannot be `TRUE` when a `data.frame` is provided as method.
#' @param fs_scores Scores ( fs) for different abundance categories of taxa associated with flow groups 1-4. A custom `data.frame` can be provided by the user.
#' @param exceptions taxa that need to be exluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' `agg` cannot be `TRUE` when a data.frame is provided as `method`.
#' @param traceB if set to `TRUE` a list as specified below will be returned.
#'
#' @keywords life
#' @importFrom dplyr '%>%' select inner_join group_by summarise rename
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble deframe
#' @importFrom stats aggregate
#'
#' @details Lotic-invertebrate Index for Flow (LIFE) was originally proposed by Extence et al. (1999). biomonitoR implements the Extence et al. (1999) version called "extence" and the version currently used in UK called "life_2017".
#' If composite is set to T the following composite families are used for extence
#'
#' \enumerate{
#'   \item Psychomyiidae (inc. Ecnomidae)
#'   \item Rhyacophilidae (inc. Glossomatidae)
#'   \item Ancylidae (inc. Acroloxidae)
#'   \item Gammaridae (inc. Crangonyctidae)
#'   \item Planariidae (inc. Dugesidae)
#'   \item Hydrobiidae (inc. Bithyniidae)
#' }
#'
#' while for "life_2017" the following are used:
#'
#'  \enumerate{
#'    \item Hydrophilidae (inc. Georissidae, Helophoridae, Hydrochidae)
#'
#'  }
#'
#' The `life` function automatically check for parent-child pairs in the scoring system, see the return section for a definition.
#' All the information used for `life` calculation can be retrieved with the function \code{\link{showscores}}.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of the `life` index.
#'  \item `taxa_df` The data.frame used for the calculation containing the abundance of taxa receiving a score.
#'  \item `abu_df` The data.frame containing scores and abundance classes for each site.
#'  \item `life_df` The data.frame used for the calculation containing scores for each site.
#'  \item `composite_taxa` Taxa aggregated following the aggregation rules when agg is not `FALSE`.
#'  \item `exceptions` A data.frame containing the containing changes made by excluding the taxa included in `exceptions`.
#'  \item `parent_child_pairs` For instance in Spanish `bmwp` both Ferrissia and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#'
#' @references Extence CA, Balbi DM, Chadd RP. 1999. River flow indexing using British benthic macroinvertebrates: a framework for setting hydroecological objectives. Regulated Rivers: Research and Management 15: 543-574.
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated LIFE scores and algorithms.

#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' life( data.agR )
#'
#' # change abundance classes
#' life( data.agR , abucl = c( 1, 9 , 99 , 999 ) )
#'
#' # provide your own score system. Scores and aggregation rules are for example purpose only.
#'
#' life_fg <- data.frame( Taxon = c( "Ephemerellidae" , "Leuctridae" , "Chironomidae" ) ,
#'  FG_Score = c( 1 , 2 , 3 ) )
#'
#' life_acc <- data.frame( Taxon = "Ephemerellidae" , Correct_Taxon = "Chironomidae" )
#'
#' fs_scores <- data.frame( FS = rep( 1:3 , each = 3 ) , ABUCLASS = rep( 1:3 , 3 ) ,
#' SCORE = c( 9 , 10 , 11 , 8 , 9 , 10 , 7 , 7 , 7 ) )
#'
#' # without aggregation rules
#' life( data.agR , method = life_fg , fs_scores = fs_scores , traceB = TRUE )
#'
#' # with aggregation
#'
#' life( data.agR , method = life_fg , agg = life_acc , fs_scores = fs_scores , traceB = TRUE )

life <- function( x , method = "extence" , abucl = c( 1 , 9 , 99 , 999 , 9999 ) , agg = FALSE , fs_scores = NULL , exceptions = NULL , traceB = FALSE ) {

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # useful for transforming data to 0-1 later
  if( inherits( x , "bin" ) ){
    BIN <- TRUE
  } else { BIN <- FALSE }

  if( ! any( identical( method , "extence" ) | identical( method , "life_2017" ) | is.data.frame( method ) ) ) stop( "method is not one of extence, life_2017 or a custom data.frame" )
  if( ! any( isFALSE( agg ) | isTRUE( agg ) | is.data.frame( agg ) ) ) stop( "agg is not one of FALSE, TRUE or a custom data.frame" )
  if( ! any( is.null( fs_scores ) | is.data.frame( fs_scores ) ) ) stop( "fs_scores is not one of NULL or a custom data.frame" )

  # Store tree for searching for inconsistencies
  Tree <- x[[ "Tree" ]][ , 1:10 ]

  numb <- c( which( names( x ) == "Tree" ) , which( names( x ) == "Taxa" ) ) # position of the Tree and Taxa data.frame in the biomonitoR object that need to be removed

  # remove Tree and Taxa data.frame
  x <- x[ -numb ]
  st.names <- names( x[[ 1 ]][ -1 ] ) # names of the sampled sites

  # create dummy variables to avoid R CMD check NOTES


  z <- w <- Taxon <- Score <- Sample <- FS <- FG_SCORE <- FG_Score <- SCORE <- NULL

  # the following if statement is to allow the users to provide their own life scores and aggregation rules.
  # y represents the method to be used

  if( is.data.frame( method ) ){
    if( ! ( isFALSE( agg ) | is.data.frame( agg ) ) ){
      stop( "When method is a data.frame agg needs to be FALSE or a data.frame containing the aggregation rules")
    }
    y <- method
    if( ! is.data.frame( fs_scores ) ) stop( "fs_scores needed when method is a data.frame provided by the user" )
    w <- fs_scores
    if( is.data.frame( agg ) ){
      z <- agg
    }

  } else {

    if( ! ( isTRUE( agg ) | isFALSE( agg ) ) ) stop( "When using the deafult method agg can only be TRUE or FALSE")
    if( ! is.null( fs_scores ) ) stop( "When using deafult methods fs_scores can only be NULL")

    # assign the default scores and aggregation rules as needed by the user
    if( identical( method , "extence" ) ){
      y <- life_scores_fam_extence
      w <- life_fs
      if( isTRUE( agg ) ){
        z <- life_acc_fam_extence
      }
    }

    if( identical( method , "life_2017" ) ){
      y <- life_scores_fam_2017
      w <- life_fs
      if( isTRUE( agg ) ){
        z <- z <- life_acc_fam_2017
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

  DF <- manage_inconsistencies( DF = DF , Tree = Tree )
  if( ! is.data.frame( DF ) ){
    incon <- DF[[ 2 ]]
    DF <- DF[[ 1 ]]
  }

  # transform the data.frame from abundance to presence-absence if needed
  if( BIN ){
    DF <- to_bin( DF )
  }

  if( traceB ){
    df2 <- DF
  }

  class.fun <- function( x ) cut( x , breaks = c( abucl , 10^18 ) , labels = 1:length( abucl ) , include.lowest = TRUE , right = TRUE )
  abu.class <- apply( apply( DF[ , - 1 , drop = FALSE ] , 2 , class.fun ) , 2 , as.numeric )
  abu.class[ is.na( abu.class ) ] <- 0
  DF <- data.frame( Taxon = DF[ , 1 ] , abu.class , check.names = FALSE )
  tot.mer <- merge( y , DF )

  res <- tot.mer %>% pivot_longer( -c( Taxon , FG_Score ) , names_to = "Sample", values_to = "ABUCLASS" ) %>%
    dplyr::rename( "FS" = "FG_Score" ) %>% inner_join( w , by = c( "FS" , "ABUCLASS" ) ) %>% select( c( "Sample" , "SCORE" ) ) %>%
    group_by( Sample ) %>% summarise( LIFE = mean( SCORE ) ) %>% deframe()

  # assure to have results if all indicator taxa are missing from a sample
  if( length( res ) != length( st.names ) ){
    temp <- merge( data.frame( Sites = st.names ) , data.frame( Sites = names( res ) , res_tot = res ) , all = TRUE )
    res <- temp$res_tot
    names( res ) <- temp$Sites
  }

  res <- res[ match( st.names , names( res ) ) ]

  if( traceB == FALSE ){
    res
  } else {
    if(  ! exists( "taxa.to.change" ) ){
      df3 <- NA
    } else { df3 <- taxa.to.change }
    if( exists( "exce" , inherits = FALSE  ) ){
      df4 <- exce
    } else { df4 <- "none" }
    if( exists( "incon" , inherits = FALSE  ) ){
      df5 <- incon
    } else { df5 <- "none" }
    temp <- tot.mer %>% pivot_longer( -c( Taxon , FG_Score ) , names_to = "Sample", values_to = "ABUCLASS" ) %>%
      dplyr::rename( "FS" = "FG_Score" ) %>% inner_join( w , by = c( "FS" , "ABUCLASS" ) ) %>% select( c( "Taxon" , "Sample" , "SCORE" ) ) %>%
      pivot_wider( names_from = Sample , values_from = SCORE , values_fill = list( SCORE = 0 ) ) %>% as.data.frame( )
    res <- list( results = res , taxa_df = df2 , abu_cl = tot.mer , life_df = temp , composite_taxa = df3 , exceptions = df4 ,  parent_child_pairs = df5 )
    res
  }
}
