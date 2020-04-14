rename <- function( x, customx = FALSE , group = "none" ){

  # use the function suggestNames to check for mispelled names and to suggest correct names
  dfName <- deparse( substitute( x ) )
  if(customx == F){
    tempNames <- suggestNames( x , group = group )
  } else {
    tempNames <- suggestNames(x, custom = T, group = group)
  }


  # if no errors were present int the user dataset return the dataset unchanged
  # otherwise change the names
  if( length( tempNames ) == 0 ){
    message("All names are correct")
    return( x )
  }  else {
  # store right and wrong names in vectors correct and wrong
  wrong <- tempNames$wrongNames
  correct <- tempNames$correctNames

  # get the taxa list from the user dataset and calculate the length of the taxa list
  result <- x$Taxa
  n <- length( wrong )

  # check if the length of correct names equals those of the wrong names
  # if not there is a problem
  if ( length( wrong ) != length( correct ) ){
    stop("pattern and replacement do not have the same length.")
  }

  # change the wrong names in the vector results that is storing the taxa list
  for(i in 1:n){
    result <- gsub( wrong[ i ] , correct[ i ] , result )
  }

  # change the taxa list of the dataset provided by the user and return the result
  x$Taxa <- result
  # the underscore is because the species names in the dictionary have the underscore
  # for example Baetis rhodani is stored as Baetis_rhodani because otherwise it would be a problem
  # to get correct species names. Probably need to be improved in the next versions of biomonitoR
  x$Taxa <- sub( "_" , " ", x$Taxa )
  x
  }
}
