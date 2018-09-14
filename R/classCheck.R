# this is a utilty function to check the class of objects provided by thr user

classCheck <- function(x, group = "none"){
  # check if the object is of class biomonitoR
  
  if( ! "biomonitoR" %in% class( x ) )(stop("The object you provided is not an object of class biomonitoR"))
  
  res <- "Maybe you have chosen the wrong taxonomic group"
  
  if( !"custom" %in% class( x ) ) {
    
    # check if the object is of class mi (macroinvertebrates)
    if(group == "mi"){
      if( ! "mi" %in% class( x ) ) stop( ( res ) )
    }
    
    # check if the object is of class mf (macrophytes)
    if(group == "mf"){
      if( ! "mf" %in% class( x ) ) stop( ( res ) )
    }
    
    # check if the object is of class fi (fish)
    if(group == "fi"){
      if( ! "fi" %in% class( x ) ) stop( ( res ) )
    }
  }
  
  
}


