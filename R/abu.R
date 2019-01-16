abu <- function(x){

  # check if the object d is of class "biomonitoR"
  classCheck(x)

  nord <- apply( x[[ "Phylum" ]][ , -1 , drop = F ], 2, FUN = sum)
  return(nord)
}

