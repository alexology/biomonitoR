#' manageTraits
#'
#' A function to select traits based on taxonomic distance.
#'
#' @param x result of the `traitScaling` function.
#' @param method can be `nearest`, `nearest+`, `nearest-`, `nearest+-` and `neareast-+`.
#'  Please see details for further information.
#' @param traceB if `TRUE` it will return a vector containing taxa excluded from the selection process
#'  because they did not meet the reuirments of the selection.
#'
#' @details Method `nearest` select the traits belonging to the nearest taxa irrispective of their position
#' in the taxonomic tree.
#' Method `nearest+` select the traits belonging to the nearest taxa that have a taxonomic level equal or finer
#' than the target one. Method `nearest-` do the opposite.
#' Method `nearest+` select the traits belonging to the nearest taxa giving priority to taxa having
#' taxonomic level equal or finer than the target one. Method `nearest-+` do the opposite.
#'
#' @export
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.ts <- traitScaling( data.agR )
#'
#' # select only the nearest traits
#' data.ts.sub <- manageTraits( data.ts , method = "nearest+-" )
#'

manageTraits <- function( x , method = "nearest+-" , traceB =  FALSE ){

  taxa <- x$Taxa
  x.split <- split( x , taxa )


  if( identical( method , "nearest" ) ) ( res <- nearest( x.split ) )
  if( identical( method , "nearest+" ) ) ( res <- nearest_pos( x.split ) )
  if( identical( method , "nearest-" ) ) ( res <- nearest_neg( x.split ) )
  if( identical( method , "nearest+-" ) ) ( res <- nearest_pos_neg( x.split ) )
  if( identical( method , "nearest-+" ) ) ( res <- nearest_neg_pos( x.split ) )


  res <- do.call( "rbind" , res )
  rownames( res ) <- NULL

  if( ! traceB ){
    res
  } else {
    missing.taxa <- taxa[ ! taxa %in% res$Taxa ]
    if( length( missing.taxa ) == 0 ){
      missing.taxa <- "none"
    }
    list( results = res , taxa_excluded = missing.taxa )
  }

}


nearest <- function( x ){
  temp.res <- lapply( x , FUN = function( z ) z[ abs( z[ , "Taxonomic_distance" ] ) >= 0 , ] )
  more_than_0 <- lapply( temp.res , FUN = function( z )( nrow( z[ abs( z[ , "Taxonomic_distance" ] ) >= 0 , ] ) )  )
  if( sum( unlist(more_than_0 ) ) == 0 ) ( stop( "Something went wrong. Please contact the mainteiner." ) )
  temp.res <-  temp.res[ more_than_0 > 0 ]
  lapply( temp.res , function( z ) z[ z[ , "Taxonomic_distance" ] == min( z[ , "Taxonomic_distance" ] )  , ] )
}

nearest_pos <- function( x ){
  temp.res <- lapply( x , FUN = function( z ) z[ z[ , "Taxonomic_distance" ] >= 0 , ]  )
  more_than_0 <- lapply( temp.res , FUN = function( z )( nrow( z[ abs( z[ , "Taxonomic_distance" ] ) >= 0 , ] ) )  )
  if( sum( unlist(more_than_0 ) ) == 0 ) ( stop( "No traits assigned at lower level") )
  temp.res <-  temp.res[ more_than_0 > 0 ]
  lapply( temp.res , function( z ) z[ z[ , "Taxonomic_distance" ] == min( z[ , "Taxonomic_distance" ] )  , ] )
}

nearest_neg <- function( x ){
  temp.res <- lapply( x , FUN = function( z ) z[ z[ , "Taxonomic_distance" ] <= 0 , ] )
  more_than_0 <- lapply( temp.res , FUN = function( z )( nrow( z[ abs( z[ , "Taxonomic_distance" ] ) >= 0 , ] ) )  )
  if( sum( unlist(more_than_0 ) ) == 0 ) ( stop( "No traits assigned at higher level") )
  temp.res <-  temp.res[ more_than_0 > 0 ]
  lapply( temp.res , function( z ) z[ z[ , "Taxonomic_distance" ] == max( z[ , "Taxonomic_distance" ] )  , ] )
}

nearest_pos_neg <- function( x ){
  temp_pos.res <- lapply( x , FUN = function( z ) nrow( z[ z[ , "Taxonomic_distance" ] >= 0 , ] ) > 0 )
  temp_neg.res <- lapply( x , FUN = function( z ) nrow( z[ z[ , "Taxonomic_distance" ] < 0 , ] ) > 0 )

  if( sum( unlist( temp_pos.res ) ) == 0 & sum( unlist( temp_neg.res ) ) == 0 ) ( stop( "Something went wrong. Please contact the mainteiner." ) )


  # assing the nearest traits, at first to positive and then to negative
  if( sum( unlist( temp_pos.res ) ) > 0 ){
    x[ unlist( temp_pos.res )  ] <- lapply( x[ unlist( temp_pos.res ) ] , FUN = function( z ) z[ z[ , "Taxonomic_distance" ] == min( z[ , "Taxonomic_distance" ] )  , ] )
  }

  if( sum( unlist( temp_neg.res ) ) > 0 ){
    x[ unlist( temp_neg.res ) ] <- lapply( x[ unlist( temp_neg.res ) ] , FUN = function( z ) z[ z[ , "Taxonomic_distance" ] == max( z[ , "Taxonomic_distance" ] )  , ] )
  }

  x

}


nearest_neg_pos <- function( x ){
  temp_neg.res <- lapply( x , FUN = function( z ) nrow( z[ z[ , "Taxonomic_distance" ] <= 0 , ] ) > 0 )
  temp_pos.res <- lapply( x , FUN = function( z ) nrow( z[ z[ , "Taxonomic_distance" ] > 0 , ] ) > 0 )


  if( sum( unlist( temp_pos.res ) ) == 0 & sum( unlist( temp_neg.res ) ) == 0 ) ( stop( "Something went wrong. Please contact the mainteiner." ) )


  # assing the nearest traits, at first to positive and then to negative
  if( sum( unlist( temp_pos.res ) ) > 0 ){
    x[ unlist( temp_pos.res )  ] <- lapply( x[ unlist( temp_pos.res ) ] , FUN = function( z ) z[ z[ , "Taxonomic_distance" ] == min( z[ , "Taxonomic_distance" ] )  , ] )
  }

  if( sum( unlist( temp_neg.res ) ) > 0 ){
    x[ unlist( temp_neg.res ) ] <- lapply( x[ unlist( temp_neg.res ) ] , FUN = function( z ) z[ z[ , "Taxonomic_distance" ] == max( z[ , "Taxonomic_distance" ] )  , ] )
  }

  x

}




