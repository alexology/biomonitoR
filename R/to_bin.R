# helper function to tranform data to presence-absence

to_bin <- function( x ){
  x_num <- unlist( lapply( x , is.numeric ) )
  x_bin <- x[ , x_num ]
  x_bin[ x_bin > 0 ] <- 1
  x.df <- data.frame( x[ , ! x_num , drop = FALSE ] , x_bin , check.names = FALSE )
  class( x.df[ , ! x_num  ] ) <- class( x[ , ! x_num   ] )
  x.df
}
