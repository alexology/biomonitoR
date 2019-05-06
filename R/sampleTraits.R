#' @describeIn traitScaling sampleTraits

sampleTraits <- function( x ){
  unique.tx <- unique( x[ , "Taxa" ] )
  res <- lapply( unique.tx ,  function( y ) ( x[ sample( rownames( x[ x$Taxa == y , ]  ) , 1 ) , ] ) )
  res <- do.call( rbind , res )
  rownames( res ) <- NULL
  res
}
