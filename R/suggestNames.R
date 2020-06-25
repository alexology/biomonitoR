#' @importFrom utils select.list
#' @importFrom hunspell dictionary hunspell_check hunspell_suggest

suggestNames <- function( x , custom = F , group = group , ... ){

  # set the path of the dictionary. If group is mi or mf the dictionary is stored in the package path
  # while if group is set to custom the dictionary is created in the user directory.
  # the check of the names relies on the package hunspell
  if( custom == F ){
    if( group == "mi" ){
      dic.path <- system.file( "dict" , "mi_dictionary.txt", package = "biomonitoR" )
      # very important to set cache equal to FALSE, otherwise suggestNames will provide inconsistent results.
      dictio <- dictionary( dic.path, cache = F )
    }
    if( group == "mf" ){
      dic.path <- system.file( "dict" , "mf_dictionary.txt" , package = "biomonitoR" )
      dictio <- dictionary( dic.path , cache = F )
    }
    if( group == "fi" ){
      dic.path <- system.file( "dict" , "fi_dictionary.txt" , package = "biomonitoR" )
      dictio <- dictionary( dic.path , cache = F )
    }
  }
  if( custom == T ){
    dic.path <- c( paste(getwd() , "/custom_dictio.dic" , sep = "" ) )
    dictio <- dictionary( dic.path , cache = F )
  }

  # get unique names present in the taxa list
  taxaCar <- unique( as.character( x$Taxa ) )


  # replace space with underscore to be compatible with hunspell
  taxaCar <- gsub( " " , '_' , taxaCar )

  # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
  # hunspell_check and hunspell_suggest are from the package hunspell
  nameCheck <- hunspell_check( taxaCar , dict = dictio )
  nameSuggest <- hunspell_suggest( taxaCar , dict = dictio )

  # This part of the function change the wrong names to correct
  # the user is provided with an interactive selection interface with select.list

  n <- which( ! nameCheck ) #number of wrong names
  if( length( n ) == 0 ){
    return( n )
  }
  else{
    temp <- rep( NA , length( n ) ) # vector to store user choices
    wrongName <- as.vector( taxaCar[ n ] ) # vector of wrong taxa names
    correctName <- nameSuggest[ n ] # list of suggestions
    for( i in 1:length( n ) ){
      choice <- c( correctName[[ i ]], "exit" ) # choices provided to the user
      # provide alternative names to the user and if the user can't find the correct name he must exit
      temp[ i ] <- select.list( choice , title = wrongName[ i ] )
      if( temp[ i ] == "exit" ) stop()
    }
    # probably it is not necessary to do Sys.unsetenv, probably dangerous if dic.path is an environmental variable itself
    # commented for now
    # Sys.unsetenv( dic.path )
    return( list( "wrongNames" = wrongName , "correctNames" = temp ) )
  }
}
