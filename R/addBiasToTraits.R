#' Add bias to duplicated fuzzy traits
#'
#' @param x results of `traitsMean` or a `data.frame` with a similar structure.
#' @param colB A vector that contains the number of modalities for each trait
#' @param SD the amount of bias.
#' @param traceB if `TRUE` returns a vector of taxa with duplicated traits.
#'
#' @details This function works by adding a small bias to duplicated fuzzy traits.
#' The bias is added with the `rnorm` function.
#'
#' @export
#' @importFrom dplyr group_by_all filter '%>%'

addBiasToTraits <- function( x , colB = NULL , SD = 0.001 , traceB = TRUE ){

  # simulate modalities of duplicated traits
  # at first a small bias is addedd to the original traits with rnorm
  # traits with 0 values are preserved by multypling biased values by ( x > 0 )
  # to be sure to have positive values we take the abolute values and data
  # are forced to sum to 1 by dividing for sum( abs( simDat ) )

  if( is.null( colB ) ) stop( "Please set colB." )

  simTraits <- function( x , ... ){
    simDat <- x + rnorm( length( x ) , ... )
    simDat <- simDat * ( x > 0 )
    abs( simDat ) / sum( abs( simDat ) )
  }


  temp.df <- x[ , ! colnames( x ) %in% "Taxa" ] %>% group_by_all %>% count %>% filter( n > 1 )

  if( nrow( temp.df ) == 0 ){
    stop( "No duplicated traits" )
  }


  # initialize vector to store row number of duplicated taxa
  names.vec <- c()

  # -ncol( temp.df ) because the last column of temp df is the duplicated count
  for( i in 1:nrow( temp.df ) ) {
    res.names  <- which( apply(  x[ , ! colnames( x ) %in% "Taxa" ] , 1 , function( z )return(all( z == t( temp.df[ i , -ncol( temp.df ) , drop = FALSE ] )[ , 1 ]   ) )  ) )
    names.vec <- c( names.vec , res.names )
  }



  names.tax <- x[ names.vec , colnames( x ) %in% "Taxa" ]

  n.rep <- temp.df$n

  res.simu <- matrix( NA, ncol = ncol( x ) - 1 , nrow = 0 )

  colnames( res.simu ) <- colnames( x )[ - 1 ]

  for( i in 1:nrow( temp.df ) ){
    # subset the ith trait to simulate
    temp.trait <- temp.df[ i , ]
    # initialize an empty matrix to store the results of the simulation of the ith trait
    empty.simu <- matrix(NA , ncol = 0 , nrow = n.rep[ i ]   )
    j2 <- 0
    for( j in 1:length( colB ) ) {
      j1 <- j2 + 1
      j2 <- j2 + colB[ j ]
      mod.simu <- temp.trait[  , j1:j2 , drop = FALSE ]
      simu.res <- t( replicate( n.rep[ i ] , simTraits( mod.simu , mean = 0 , sd = SD ) , simplify = TRUE  ) )
      empty.simu <- cbind( empty.simu , simu.res )
    }
    res.simu <- rbind( res.simu , empty.simu )
  }

  res.simu <- as.data.frame( sapply( as.data.frame( res.simu ), as.numeric ) )
  res.simu <- data.frame( Taxa = names.tax , res.simu )
  res.simu <- rbind( x[ -names.vec , ] , res.simu )
  res.simu <- res.simu[ match( x[ , colnames( x ) %in% "Taxa" ] , res.simu[ , colnames( res.simu ) %in% "Taxa" ]  ) , ]
  rownames( res.simu ) <- NULL

  if( ! traceB ){
    res.simu
  } else {
    list( results = res.simu , duplicated_traits = as.character( names.tax )  )
  }


}








