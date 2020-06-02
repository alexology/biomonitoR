# this is a utilty function to check the class of objects provided by the user

classCheck <- function( x , group = "none" ){
  # check if the object is of class biomonitoR

  if( ! "biomonitoR" %in% class( x ) )( stop( "The object you provided is not an object of class biomonitoR" ) )

  res <- "It seems that you used your own reference database. Please check the consistency of the taxonomy used for calculating the index with those of your reference database to have reliable results."

  if( !"custom" %in% class( x ) ) {

    # check if the object is of class mi (macroinvertebrates)
    if( group == "mi" ){
      if( ! "mi" %in% class( x ) ) warning( ( res ) )
    }

    # check if the object is of class mf (macrophytes)
    if( group == "mf" ){
      if( ! "mf" %in% class( x ) ) warning( ( res ) )
    }

    # check if the object is of class fi (fish)
    if( group == "fi" ){
      if( ! "fi" %in% class( x ) ) warning( ( res ) )
    }
  }
}


