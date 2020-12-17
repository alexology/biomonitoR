#' Functional dispersion
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates the functional richness based on trait categories.
#'
#'
#'
#' Functional richness (FRic) represents the amount of functional space filled by
#' the community (Villeger et al., 2008) and it is related to the community use of
#' resources and productivity (Mason et al., 2005). FRic is defined by the trait
#' extremes and thus reflects the potential maximum functional dissimilarity.
#' FRic is calculated as the hypervolume enclosing the functional space filled
#' by the community. For this the convex hull approach is applied (Cornwell et al.,
#' 2006) using the Quickhull algorithm (Barber et al., 1996) that estimates the minimum
#'  convex hull which includes all the species considered in the previously defined
#'  functional space. Basically, this algorithm determines the taxa in the most extreme
#'  points of the functional space, links them to build the convex hull in order to
#'  calculate the volume inside it. In particular, the convex hull of a set of points
#'  S in n dimensions is the intersection of all convex sets containing S.
#'  For N points , ..., , the convex hull C is then given by the expression:
#'
#'  \deqn{C = \sum_{j=1}^{N} \lambda_j \ p_j  : \ \lambda_j \geq \ for \ all \ j \ and \ \sum_{j=1}^{N} \lambda_j = 1}
#'
#' The functional T dimensional space is built using a certain number of dimensions
#' (T) determined by the axes of a principal component analysis based on the trait
#'  dissimilarity matrix. Using a poor quality functional space could led to
#'  a biased assessment of FRic and false ecological conclusions so the number of
#'  axes retained for its estimation is case-specific and decided following the
#'  method proposed in Maire et al., (2015): a pragmatic approach consisting of
#'  computing all the possible functional spaces and selecting the most
#'  parsimonious one. This method uses the mSD index (mean squared deviation
#'  between the initial functional distance and the scaled distance in the
#'  functional space), which accounts explicitly for the deviation between
#'  the initial and final distance and penalizes the strong deviation. The mSD
#'  index has been widely used in statistics to assess errors and it has been
#'  demonstrated to work in different contexts and situations (Maire et al., 2015).
#'  In addition, when using Gower's distance the mSD ranges from 0 and 1,
#'  which helps to interpret quality. Finally, the resulting FRic variable
#'  is standardized by its maximum, ranging from 0 to 1. In addition, it must
#'  be considered that the number of taxa must be higher than the number of traits
#'  to have reliable FRic values (Villeger et al., 2008).
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
#' @param nbdim number of dimensions for the multidimensional functional spaces.
#' We suggest to keep `nbdim` as low as possible.
#' By default `biomonitoR` set the number of dimensions to 2. Select `auto` if you want the automated selection
#' approach according to Maire et al. (2015).
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
#' @return a vector with fuzzy functional richness results.
#' \enumerate{
#'  \item **results**: results of the ffred function;
#'  \item **traits**: a data.frame containing the traits used for the calculations;
#'  \item **taxa**: a data.frame conaining the taxa used for th calculations;
#'  \item **nbdim**: number of dimensions used after calculatin the quality of functional spaces according to Maire et al., (2015);
#'  \item **correction**: the type of correction used.
#'  \item **NA_detection**: a data.frame containing taxa on the first column and the corresponding trais with NAs on the second column.
#'  \item **duplicated_traits**: if present, list the taxa with the same traits.
#' }
#'
#' @importFrom ade4 prep.fuzzy dudi.pco is.euclid cailliez lingoes prep.binary prep.circular bicenter.wt
#' @importFrom graphics abline barplot par
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
#' f_disp( data.agR , traitDB = data.ts.av , type = "F" , colB = colB )
#' f_disp( data.agR , traitDB = data.ts.av , type = "F" , colB = colB ,
#'        nbdim = 10 , correction = "cailliez" )
#'
#' library( ade4 )
#'
#' rownames( data.ts.av ) <- data.ts.av$Taxa
#' traits.prep <- prep.fuzzy( data.ts.av[ , -1 ], col.blocks = colB )
#'
#' traits.dist <- ktab.list.df( list( traits.prep ) )
#' traits.dist <- dist.ktab( traits.dist , type = "F" )
#'
#' f_disp( data.agR , traitDB = traits.dist   )
#'
#'
#' @seealso [aggregatoR]
#'
#' @references Barber, C. B., Dobkin, D. P., & Huhdanpaa, H. (1996).
#'   The quickhull algorithm for convex hulls. ACM Transactions on
#'   Mathematical Software (TOMS), 22(4), 469-483.
#' @references Cornwell, W. K., Schwilk, D. W., & Ackerly, D. D. (2006).
#'   A trait-based test for habitat filtering: convex hull volume.
#'   Ecology, 87(6), 1465-1471
#' @references Maire, E., Grenouillet, G., Brosse, S., & Villeger, S. (2015).
#'   How many dimensions are needed to accurately assess functional diversity?
#'   A pragmatic approach for assessing the quality of functional spaces. Global
#'   Ecology and Biogeography, 24(6), 728-740.
#' @references Mason, N. W., Mouillot, D., Lee, W. G., and
#'   Wilson, J. B. (2005). Functional richness, functional evenness and functional
#'   divergence: the primary components of functional diversity. Oikos, 111(1),
#'   112-118.
#' @references Villeger, S., Mason, N. W., & Mouillot, D.
#'   (2008). New multidimensional functional diversity indices for a
#'   multifaceted framework in functional ecology. Ecology, 89(8), 2290-2301.
#'
#' @export

f_disp <- function( x , traitDB = NULL, taxLev = "Taxa" , type = NULL , traitSel = FALSE , colB = NULL,  nbdim = 2 , distance = "gower", zerodist_rm = FALSE , correction = "none" , traceB = FALSE , set_param = list( max_nbdim = 15 , prec = "Qt" , tol = 1e-07 , cor.zero = TRUE ) ){

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
    MES <- "At least a pair of species has the same traits. Depending on your needs, this could be an issue."
    message( MES )
  } else {
    MES <- "no taxa with the same traits"
  }

  if( ! euclid.dist.mat )( message( "Negative eigenvalues found, please consider to add a correction to the pcoa" ) )


  if( identical( nbdim , "auto" ) ){
    qual_fs <- qfs( mat_dissim , nbdim = set_param$max_nbdim  )
    m.check <- qual_fs$meanSD < 0.01

    if( ! any( na.omit( m.check ) ) ) { stop ("there is no optimal number of dimension, please increase the number of dimensions" )}
    m <- min( which( m.check ) ) + 1
  } else{
    m <- nbdim
  }


  # transpose DF to fit with the fric_3d function

  abu.t <- t( DF[ , -1 ] )
  colnames( abu.t ) <- DF$Taxon

  res <- fdisp_k( mat_dissim , abu.t , m )

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
    } else { df1 <- MES }


    rownames( DF ) <- NULL

    res.list <- list( res , traitDB , DF , m , correction = correction, tax.na , df1  )
    names( res.list ) <- c( "results" , "traits" , "taxa" , "nbdim" , "correction" , "NA_detection" , "duplicated_traits"  )
    return( res.list )
  }

}
