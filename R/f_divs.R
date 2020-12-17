#' Functional diversity
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates the functional diversity.
#'
#' Rao quadratic entropy (Q; Rao, 1982) is used to estimate Functional Diversity
#' because it has been considered more appropriate than other indices (Botta-Dukat,
#' 2005; Ricotta, 2005). In this formula (Q) dij is the dissimilarity (ranging
#' from 0 to 1), between species i and j based on a set of specified functional
#' traits (i.e. effect traits, see below). This index is standardized by the maximum value to constrain the values
#' within the range of 0-1. Rao index is estimated using presence or abundance data
#' and the Euclidean transformed version of the traits-based Gower dissimilarity
#' matrix. For this, Gower's dissimilarity index (which ranges from 0 to 1) is used
#' because it can deal with traits of different nature and measuring scales
#' (continuous, nominal, binary, ordinal, etc.; see Podani 1999 for more information)
#'
#' \deqn{Q = \sum_{i=1}^{S} \sum_{j=1}^{S} d_{ij} \ p_i \ p_j}
#'
#'
#' @param x results of function aggregatoR
#' @param traitDB a trait database. Can be a `data.frame` ot a `dist` object.
#' Taxonomic level of the traits database must match those of the taxonomic database.
#' No automatic check is done by the `function`.
#' @param taxLev character string giving the taxonomic level used to retrieve
#' trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#' `"Family"` as returned by the [aggregatoR] function.
#' @param type the type of variables speciefied in `traitDB`.
#' Must be one of `F`, fuzzy, or `C`, continuous.
#' If more control is needed please consider to provide `traitDB` as a `dist` object.
#' It works only when `traitDB` is a `data.frame`, otherwise ingored.
#' @param traitSel interactively select traits.
#' @param colB A vector that contains the number of modalities for each trait.
#' Not needed when `euclidean` distance is used.
#' @param distance to be used to compute functional distances, `euclidean` or `gower`. Default to `gower`.
#' @param zerodist_rm If `TRUE` aggregates taxa with the same traits.
#' @param correction Correction methods for negative eigenvalues, can be one of `none`, `lingoes`, `cailliez`, `sqrt` and `quasi`.
#' Ignored when type is set to `C`.
#' @param traceB if `TRUE` ffrich will return a list as specified in details.
#' @param set_param a list of parameters for fine tuning the calculations.
#' `max_nbdim` set the maximum number of dimension for evaluating the quality of the functional space.
#' `prec` can be `Qt` or `QJ`, please refere to the `convhulln` documentation for more information.
#' Deafault to `QJ`, less accurate but less prone to errors.
#' `tol` a tolerance threshold for zero, see the function `is.euclid`, `lingoes` and `cailliez` from the `ade4` for more details. Default to 1e-07.
#' `cor.zero` = `TRUE` if TRUE, zero distances are not modified. see the function `is.euclid`, `lingoes` and `cailliez` from the `ade4` for more details. Default to `TRUE`.
#'
#' @details Taxa without traits assigned in the trait database are removed from both the trait and abundance databases.
#'
#' @return a vector with fuzzy functional richness results.
#' \enumerate{
#'  \item **results**: results of the ffred function;
#'  \item **traits**: a data.frame containing the traits used for the calculations;
#'  \item **taxa**: a data.frame conaining the taxa used for th calculations;
#'  \item **correction**: the type of correction used.
#'  \item **NA_detection**: a data.frame containing taxa on the first column and the corresponding trais with NAs on the second column.
#'  \item **duplicated_traits**: if present, list the taxa with the same traits.
#' }
#'
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#' @importFrom stats complete.cases na.omit
#' @importFrom ade4 ktab.list.df dist.ktab prep.fuzzy divc quasieuclid is.euclid
#' @importFrom stats dist weighted.mean
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.ts <- traitScaling( data.agR )
#' # averaging
#' data.ts.av <- traitsMean( data.ts )
#'
#' colB <- c( 8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8 )
#'
#' f_divs( data.agR , traitDB = data.ts.av , type = "F" , colB = colB )
#' f_divs( data.agR , traitDB = data.ts.av , type = "F" , colB = colB ,
#'        correction = "cailliez" )
#'
#' library( ade4 )
#'
#' rownames( data.ts.av ) <- data.ts.av$Taxa
#' traits.prep <- prep.fuzzy( data.ts.av[ , -1 ], col.blocks = colB )
#'
#' traits.dist <- ktab.list.df( list( traits.prep ) )
#' traits.dist <- dist.ktab( traits.dist , type = "F" )
#'
#' f_divs( data.agR , traitDB = traits.dist   )
#'
#'
#'
#' @seealso [aggregatoR]
#'
#' @references Biggs, R., Schluter, M., Biggs, D., Bohensky, E. L., BurnSilver, S.,
#'   Cundill, G., ... & Leitch, A. M. (2012). Toward principles for enhancing the
#'   resilience of ecosystem services. Annual Review of Environment and Resources,
#'   37, 421-448.
#' @references Botta-Dukat, Z. (2005). Rao's quadratic entropy as a measure of
#'   functional diversity based on multiple traits. Journal of Vegetation Science,
#'   16(5), 533-540.
#' @references de Bello, F., Leps, J., Lavorel, S., & Moretti, M. (2007).
#'   Importance of species abundance for assessment of trait composition:
#'   an example based on pollinator communities. Community Ecology, 8(2), 163-170.
#' @references Elmqvist, T., Folke, C., Nystrom, M., Peterson, G.,
#'   Bengtsson, J., Walker, B., & Norberg, J. (2003). Response diversity, ecosystem
#'   change, and resilience. Frontiers in Ecology and the Environment, 1(9), 488-494.
#' @references Guillemot, N., Kulbicki, M., Chabanet, P., & Vigliola, L. (2011).
#'   Functional redundancy patterns reveal non-random assembly rules in a
#'  species-rich marine assemblage. PLoS One, 6(10), e26735.
#' @references Hevia, V., Martin-Lopez, B., Palomo, S., Garcia-Llorente, M., de Bello, F.,
#'   & Gonzalez, J. A. (2017). Trait-based approaches to analyze links between the drivers
#'   of change and ecosystem services: Synthesizing existing evidence and future challenges.
#'   Ecology and evolution, 7(3), 831-844.
#' @references Hooper, D. U., Chapin, F. S., Ewel, J. J., Hector, A., Inchausti,
#'   P., Lavorel, S., et al. (2005). Effects of biodiversity on ecosystem
#'   functioning: a consensus of current knowledge. Ecological Monographs,
#'   75(1), 3-35.
#' @references Lawton, J.H. & Brown, V.K. (1993) Redundancy in ecosystems.
#'   Biodiversity and Ecosystem Function (eds E.-D. Schulze & H.A. Mooney),
#'   pp. 255-270. Springer-Verlag, Berlin.
#' @references Pillar, V. D., Blanco, C. C., Muller, S. C., Sosinski, E. E.,
#'   Joner, F., & Duarte, L. D. (2013). Functional redundancy and stability
#'   in plant communities. Journal of Vegetation Science, 24(5), 963-974.
#' @references Podani, J. (1999). Extending Gower's general coefficient of
#'   similarity to ordinal characters. Taxon, 331-340.
#' @references Rao, C. R. (1982). Diversity and dissimilarity coefficients:
#'   a unified approach. Theoretical population biology, 21(1), 24-43.
#' @references Ricotta, C. (2005). A note on functional diversity measures.
#'   Basic and Applied Ecology, 6(5), 479-486.
#' @references Rosenfeld, J. S. (2002). Functional redundancy in ecology
#'   and conservation. Oikos, 98(1), 156-162.
#' @references Schmera, D., Heino, J., Podani, J., Eros, T., & Doledec, S. (2017).
#'   Functional diversity: a review of methodology and current knowledge in
#'   freshwater macroinvertebrate research. Hydrobiologia, 787(1), 27-44.
#' @references Walker, B. H. (1992). Biodiversity and ecological redundancy.
#'   Conservation biology, 6(1), 18-23.
#'
#' @export

f_divs <- function( x , traitDB = NULL, taxLev = "Taxa" , type = NULL , traitSel = FALSE , colB = NULL,  distance = "gower", zerodist_rm = FALSE , traceB = FALSE , correction = "none" ,  set_param = NULL ){

  #  check if the object x is of class "biomonitoR"
  classCheck( x )

  # useful for transforming data to 0-1 later
  if( inherits( x , "bin" ) ){
    BIN <- TRUE
  } else { BIN <- FALSE }


  # list set_param for fine tuning of some functions
  if( is.null( set_param ) ){
    set_param <- list( max_nbdim = 15 , prec = "Qt" , tol = 1e-07 , cor.zero = TRUE )
  } else {
    set_param_def <- list( max_nbdim = 15 , prec = "Qt" , tol = 1e-07 , cor.zero = TRUE )
    set_param_def[ names( set_param ) ] <- set_param
    set_param <- set_param_def
  }

  if( is.null( traitDB ) ) stop( "Please provide traitDB" )

  if( ! is.data.frame( traitDB ) & ! class( traitDB ) %in% "dist"  ) stop( "traitDB must be a data.frame or a dist object" )

  if( is.null( type ) & is.data.frame( traitDB ) ) stop( "Please specify a type when traitDB is a data.frame" )

  if( ! identical( type , "F" ) & ! identical( type , "C" ) & is.data.frame( traitDB )  ) stop( "type must be C or F when traitDB is a data.frame" )

  if( identical( type , "C" ) & identical( distance , "gower" ) ) ( warning( "Are you sure to use gower distance when type is C?" ) )

  if( identical( type , "F" ) & identical( distance , "euclidean" ) ) ( warning( "Are you sure to use euclidean distance when type is F?" ) )


  if( is.data.frame( traitDB ) ){

    # trim and capitalise the first letter in the DB provided by the user
    traitDB[ , "Taxa"] <- as.factor( sapply( trim( traitDB[ , "Taxa"] ), capWords, USE.NAMES = F ) )

    if( identical( type , "F" ) ){
      if( is.null( colB ) ) ( stop( "Please provide colB" ) )
      # check if the number of traits in traitDB equals the sum of colB, otherwise stop
      if( ( ncol( traitDB ) - 1 ) != sum( colB ) ) ( stop( "The number of traits in traitDB is not equal to the sum of colB" ) )
    }
  }


  if( traitSel & is.data.frame( traitDB ) ){
    Index <- rep( 1:length( colB ) , colB )
    rma <- select.list( names( traitDB[ -which( names( traitDB ) %in% "Taxa")] ) , title = "Traits selection"  , graphics = TRUE , multiple = T )
    # new colB based on user trait selection, -1 because there is the column called Taxa
    colB <- as.vector( table( Index[ which( names( traitDB ) %in% rma ) - 1 ] ) )

    #  trait must have at least two modalities
    if( any( colB < 2 ) ) ( stop( "a trait must have at least two modalities" ) )

    traitDB <- traitDB %>%
      select( c("Taxa", rma) )
    # trim and capitalise the column Taxa of the user' trait database
    traitDB$Taxa <- apply( as.data.frame( trim( traitDB$Taxa ) ) , 1 , capWords)

  }



  st.names <- names( x[[ 1 ]][ -1 ] )

  DF <- x[[ taxLev ]]
  names( DF )[ 1 ] <- "Taxon"

  taxa <- as.character( DF$Taxon )
  DF$Taxon <- as.character( DF$Taxon )


  if( is.data.frame( traitDB ) ) {

    traitDB$Taxa <- as.character( traitDB$Taxa )
    names( traitDB )[ names( traitDB ) %in% "Taxa" ] <- "Taxon"

    # be sure that taxonomic and functional database have the same order and taxa
    DF <- merge( DF , traitDB[ , "Taxon" , drop = FALSE ] , by = "Taxon" )


    # transform the data.frame from abundance to presence-absence if needed
    if( BIN ){
      DF <- to_bin( DF )
    }

    traitDB <- merge( traitDB , DF[ , "Taxon" , drop = FALSE ] , by = "Taxon" )

    if( any( ! DF$Taxon == traitDB$Taxon ) ) stop( "Taxonomic and traits taxa does not match, ask the maintainer" )

    # just to be sure we are doing the right things
    rownames( traitDB ) <- traitDB$Taxon

    if( identical( type , "F" ) ) ( tr_prep <- prep.fuzzy( traitDB[ , -1 ], col.blocks = colB ) )
    if( identical( type , "B" ) ) ( tr_prep <- prep.binary( traitDB[ , -1 ], col.blocks = colB ) )
    if( identical( type , "C" ) ) ( tr_prep <- traitDB[ , -1 ] )

    rownames( tr_prep ) <- traitDB$Taxon


    # computing functional dissimilarity between species given their traits values
    if ( identical( distance  ,"gower" ) ) {
      mat_dissim <- ktab.list.df( list( tr_prep ) )
      mat_dissim <- dist.ktab( mat_dissim , type = "F")
    }

    if ( identical( distance  , "euclidean" ) ) mat_dissim <- dist( scale( tr_prep ) ) # scaling if continuous traits
  }

  if( class( traitDB ) %in% "dist" ){

    # keep only common taxa between taxonomic and traits database
    # to do this the dist object is transformed into matrix
    traitDB <- as.matrix( traitDB )
    DF <- merge( DF , data.frame( Taxon = rownames( traitDB ) ) , by = "Taxon" )


    # transform the data.frame from abundance to presence-absence if needed
    if( BIN ){
      DF <- to_bin( DF )
    }

    traitDB <- traitDB[ rownames( traitDB ) %in% DF$Taxon ,  ]
    traitDB <- traitDB[ , colnames( traitDB ) %in% DF$Taxon , drop = FALSE  ]
    traitDB <- traitDB[ match( DF$Taxon , rownames( traitDB ) ) , ]
    traitDB <- traitDB[ , match( DF$Taxon , colnames( traitDB ) ) , drop = FALSE ]

    # check if names are in the same order, both on rows and columns
    if( any( ! DF$Taxon == rownames( traitDB ) ) ) stop( "Taxonomic and traits taxa does not match, ask the maintainer" )
    if( any( ! DF$Taxon == colnames( traitDB ) ) ) stop( "Taxonomic and traits taxa does not match, ask the maintainer" )

    mat_dissim <- as.dist( traitDB )


  }


  if( zerodist_rm & any( mat_dissim < set_param$tol ) ){
    zero_corr <- zero_dist_traits( x = DF , mat_dissim = mat_dissim , BIN = BIN  )
    DF <- zero_corr[[ 1 ]]
    mat_dissim <- zero_corr[[ 2 ]]
    df1 <- zero_corr[[3]]
  }


  if( identical( distance , "gower" ) ){
    if( identical( correction , "cailliez" ) ) mat_dissim <- suppressWarnings( cailliez( mat_dissim , tol = set_param$tol , cor.zero = set_param$cor.zero ) )
    if( identical( correction , "lingoes" ) ) mat_dissim <- suppressWarnings( lingoes( mat_dissim  , tol = set_param$tol , cor.zero = set_param$cor.zero ) )
    if( identical( correction , "sqrt" ) ) mat_dissim <- suppressWarnings( sqrt( mat_dissim ) )
    if( identical( correction , "quasi" ) ) mat_dissim <- suppressWarnings( quasieuclid( mat_dissim ) )
  }

  suppressWarnings( euclid.dist.mat <- is.euclid( mat_dissim , tol = set_param$tol ) )

  if( any( mat_dissim < set_param$tol ) ){
    message( "At least a pair of species has the same traits. Depending on your needs, this could be an issue." )
  }

  # transpose DF to fit with the fric_3d function

  rownames( DF ) <- DF[ , "Taxon" ]

  suppressWarnings( raoQ <- divc( DF[ , -1 ]  , mat_dissim , scale = T)$diversity )

  if( ! euclid.dist.mat ){
    warning( "Non euclidean trait distance. Euclidean property is expected. Please use the correction options
             otherwise consider to remove taxa with the same traits." )
  }

  res <- raoQ
  names( res ) <- st.names

  if( ! traceB ){
    return( res )
  }

  if( traceB ){
    # chech for NA, it could happen that a trait is filled with NAs
    # but this can be done only when traitDB is a data.frame
    if( is.data.frame( traitDB ) ){
      if( any( is.na( tr_prep ) ) ){
        tax.na <- as.data.frame( which( is.na( tr_prep ) , arr.ind = TRUE ) )
        tax.na[ , 1 ] <- traitDB[ tax.na[ , 1 ] , 1 ]
        tax.na[ , 2 ] <- colnames( traitDB[ , -1 ] )[ tax.na[ , 2 ]  ]
        colnames( tax.na ) <- c( "Taxa" , "Traits" )
        tax.na <- tax.na[ order( tax.na[ , 1 ] ) ,  ]
        rownames( tax.na ) <- NULL
      } else { tax.na <- "No NAs detected" }

    } else {
      tax.na <- "NAs cannot be detected when traitDB is a dist object"
    }

    # prepare traits to be returned
    if( ! is.data.frame( traitDB ) ){
      # returns the distance matrix used for the calculation as a dist object
      traitDB <- as.dist( traitDB )
    }

    # prepare traits to be returned
    if( is.data.frame( traitDB ) ){
      # returns the distance matrix used for the calculation as a dist object
      rownames( traitDB ) <- NULL
    }

    if(  exists( "df1" , inherits = FALSE  ) ){
        df1 <- df1
    } else { df1 <- "no taxa with the same traits"}


    rownames( DF ) <- NULL

    res.list <- list( res , traitDB , DF ,  correction = correction , tax.na , df1 )
    names( res.list ) <- c( "results" , "traits" , "taxa" , "correction" , "NA_detection" , "duplicated_traits" )
    return( res.list )
  }
}
