#' whpt
#'
#' This function calculates *Whalley Hawkes Paisley Trigg* index according to version used in UK in 2017.
#' @param x result of the function `aggregatoR`.
#' @param method the only choice is `uk`.Users can provide their own data.frame (see examples) with a column called *Taxon*, a column called *ABUCLASS* and a column called *Scores*.
#' @param type presence only `po` or abundance `ab`.
#' @param metric possible choices are `aspt`, `ntaxa`, `bmwp`.
#' @param agg this option allows the composite family approach. It can be `FALSE`, `TRUE` or a `data.frame`.
#' If `FALSE` no aggregation will be performed, while if `TRUE` aggregation will be performed according to the rules described in Details.
#' A `data.frame` containing the aggregation rules can be provided by the user.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' `agg` cannot be `TRUE` when a `data.frame` is provided as method.
#' @param abucl Log abundance categories. Treshold are set to 1, 9, 99 and 999.
#' @param exceptions taxa that need to be exluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' `agg` cannot be `TRUE` when a data.frame is provided as `method`.
#' @param traceB if set to `TRUE` a list as specified below will be returned.
#' @keywords whpt
#' @details WHPT is a revision of BMWP and it takes into account the abundances of organisms. The following aggregation is used if agg is set equal to `TRUE`:
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
#' The `whpt` function automatically check for parent-child pairs in the scoring system, see the return section for a definition.
#' All the information used for `whpt` calculation can be retrieved with the function \code{\link{showscores} }.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of the `whpt` index.
#'  \item `taxa_df` The data.frame used for the calculation containing the abundance of taxa receiving a score.
#'  \item `abu_df` The data.frame containing fssr scores and abundance classes for each site.
#'  \item `whpt_df` The data.frame used for the calculation containing scores for each site.
#'  \item `composite_taxa` Taxa aggregated following the aggregation rules when agg is not `NULL`.
#'  \item `exceptions` A data.frame containing the con.
#'  \item `parent_child_pairs` For instance in Spanish `aspt` both Ferrissia and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated WHPT scores and algorithms.
#' @importFrom dplyr '%>%' select inner_join group_by summarise filter
#' @importFrom tidyr pivot_longer
#' @importFrom tibble deframe
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{asBiomonitor}}, \code{\link{aspt}}, \code{\link{bmwp}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' whpt( data.agR )
#' whpt( data.agR , metric = "bmwp" )
#' whpt( data.agR , type = "po"  , metric = "bmwp" )
#'
#' # take a look to the metrics used for whpt calculation
#' # only the first 6 rows of each database are shown
#'
#' lapply( showscores( "whpt" , "uk" ) , head )

whpt <- function(x, method = "uk" , type = "ab", metric = "aspt", agg = FALSE , abucl = c( 1 , 9 , 99 , 999 ) , exceptions = NULL , traceB = FALSE  ){

  # check if the object x is of class "biomonitoR"
  classCheck( x )

  # useful for transforming data to 0-1 later
  if( inherits( x , "bin" ) ){
    BIN <- TRUE
  } else { BIN <- FALSE }


  if( ! identical( method , "uk" ) & !is.data.frame( method )  ){
    stop( "Method need to be set to uk or a custom data.frame" )
  }

  if( ! identical( type , "po" ) & ! identical( type , "ab" ) ){
    stop( "Please provide a valide type: po or ab" )
  }

  if( ! identical( metric , "aspt" ) & ! identical( metric , "ntaxa" )  & ! identical( metric , "bmwp" ) ){
    stop("Please provide a valide metric: aspt, ntaxa or bmwp")
  }


  # Store tree for searching for inconsistencies
  Tree <- x[[ "Tree" ]][ , 1:10 ]

  numb <- c( which( names( x ) == "Tree" ) , which( names( x ) == "Taxa" ) ) # position of the Tree and Taxa data.frame in the biomonitoR object that need to be removed

  # remove Tree and Taxa data.frame
  x <- x[ -numb ]
  st.names <- names( x[[ 1 ]][ -1 ] ) # names of the sampled sites

  # create dummy variables to avoid R CMD check NOTES


  z <- Abundance <- Taxon <- Sample <- Scores <- ABUCLASS <- NULL


  # the following if statement is to allow the users to provide their own psi scores and aggregation rules.
  # y represents the method to be used

  # the following if statement is to allow the users to provide their own bmwp scores and aggregation rules.
  # y represents the method to be used


  if( is.data.frame( method ) == TRUE ){
    if( ! ( isFALSE( agg ) | is.data.frame( agg ) ) ){
      stop( "When method is a data.frame agg needs to be FALSE or a data.frame containing the aggregation rules")
    }

    if( isFALSE( agg ) ){
      y <- method
    } else {
      y <- method
      z <- agg
    }
  } else {

    if( ! ( isTRUE( agg ) | isFALSE( agg ) ) ) stop( "When using the deafult method agg can only be TRUE or FALSE")

    # assign the default scores and aggregation rules as needed by the user

    if( identical( method , "uk" ) ) {
      y <- whpt_scores_fam_uk

      if( isTRUE( agg ) ){
        z <- whpt_acc_fam_uk
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

  # WHPT scores have a long format, unsuitable for merging with DF
  # make a temp data.frame to store unique taxa names.

  y.temp <- data.frame( Taxon = unique( y[ , "Taxon" ] )  )

  # merge the new data.frame with the score data.frame and change
  # the names of the taxa according to the aggregation rules if needed
  DF <- merge( DF, y.temp[ , "Taxon" , drop = FALSE ] )
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
  abu.class <- apply( apply( DF[ , -1 , drop = FALSE ] , 2 , class.fun ) , 2 , as.numeric )
  abu.class[ is.na( abu.class ) ] <- 0
  tot.mer <- data.frame( DF[ , 1 , drop = FALSE ] , abu.class , check.names = FALSE )

  # to avoid warning from inner_join
  y$Taxon <- as.character( y$Taxon )
  tot.mer$Taxon <- as.character( tot.mer$Taxon )

  if( identical( type , "ab" ) & BIN ) stop( "Type cannot be set to abundance if presence-absence data are used")

  if( type == "po" ){
    y <- subset( y , ABUCLASS == -1 )
    res.tot <- tot.mer %>% pivot_longer( -c( Taxon ) , names_to = "Sample", values_to = "Abundance" ) %>%
      filter( Abundance > 0 ) %>% inner_join( y , by = c( "Taxon" ) ) %>% select( c( "Sample" , "Scores" ) )
    whpt_df <- tot.mer %>% pivot_longer( -c( Taxon ) , names_to = "Sample", values_to = "Abundance" ) %>%
      filter( Abundance > 0 ) %>% inner_join( y , by = c( "Taxon" ) ) %>% as.data.frame( )
  }

  if(type == "ab"){
    y <- subset( y , ABUCLASS != -1 )
    res.tot <- tot.mer %>% pivot_longer( -c( Taxon ) , names_to = "Sample", values_to = "ABUCLASS" ) %>%
      inner_join( y , by = c( "Taxon" , "ABUCLASS" ) )
    whpt.df <- res.tot %>% as.data.frame( )
  }

  if(metric == "aspt"){
    res <- res.tot %>% group_by( Sample ) %>% summarise( mean( Scores ) ) %>% deframe()
  }
  if(metric == "ntaxa"){
    res <- res.tot %>% group_by( Sample ) %>% summarise( length( Scores ) ) %>% deframe()
  }
  if(metric == "bmwp"){
    res <- res.tot %>% group_by( Sample ) %>% summarise( sum( Scores ) ) %>% deframe()
  }

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
    if( ! exists( "taxa.to.change" , inherits = FALSE ) ){
      df3 <- "none"
    } else { df3 <- taxa.to.change }
    if( exists( "exce" , inherits = FALSE  ) ){
      df1 <- exce
    } else { df4 <- "none" }
    if( exists( "incon" , inherits = FALSE  ) ){
      df5 <- incon
    } else { df5 <- "none" }

    res <- list( results = res , taxa_df = df2 , abu_cl = tot.mer , whpt_df = whpt_df , composite_taxa = df3 ,  exceptions = df4 ,  parent_child_pairs = df5 )
    res
  }

}
