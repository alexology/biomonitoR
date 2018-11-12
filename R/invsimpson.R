#' @describeIn allindices Inverse Simpson

invsimpson <- function(x, taxLev = "Family"){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  df <-  x[[taxLev]]
  if("unassigned" %in% df[ , 1]){
    z <- which(df[ ,1 ] == "unassigned")
    df<- df[ -z ,] # remove unassigned row from the species count
  }

  res <- apply(df[ , -1, drop = FALSE ], 2, FUN = function(x){Pi( x, index = "Invsimpson")})

  return( res )
}
