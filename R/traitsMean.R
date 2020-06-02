#' @describeIn traitScaling average traits values for fuzzy data
#' @importFrom ade4 prep.fuzzy.var

traitsMean <- function( x , colB = NULL  ){

  if( is.null( colB ) ) {
    colB <- c( 8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8 )
  } else {
    colB <- colB
  }

  x <- x[ , -c( 2:5 ) , drop = FALSE ]
  x[ is.na( x ) ] <- 0
  trait_mean <- aggregate( . ~ Taxa , data = x , FUN = mean  )
  rownames( trait_mean ) <- trait_mean[ , 1 ]
  trait_mean <- trait_mean[ , - 1 , drop = FALSE ]
  trait_res <- prep.fuzzy( trait_mean , col.blocks = colB )
  rownames( trait_res ) <- NULL
  res <- data.frame( Taxa = as.factor( rownames( trait_mean ) ) , trait_res )
  res
}


