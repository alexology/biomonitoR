#' quickRename
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("deprecated") }
#'
#' This function allow the user to change the wrong taxa names according to the deault reference databases implemented in `biomonitoR`.
#'
#' @param x a data.frame as specified in details
#' @param group biotic group of interest. Possible values are `mi` for macroinvertebrates and `mf` for macrophytes. Default to `mi`.
#' @param write.table if `TRUE` `quickRename` will save a txt file with changes provided by the user
#'
#' @keywords asBiomonitor
#'
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' The function `quickRename` will suggest correct name while allowing the user to insert a name (Enter taxon name).
#' @export
#' @seealso \code{\link{asBiomonitor}}


quickRename <- function( x , group = "mi", write.table = FALSE ){

  .Deprecated( "quick_traits" , package = "biomonitoR" )

  if( is.null( group ) == TRUE ){
    stop("Please provide a valid group")
  }

  if( ! any( names( x ) %in% "Taxa" ) ) ( stop( " A column called Taxa is missing" ) )

  dfName <- deparse( substitute( x ) )
  tempNames <- suggestUserNames( x , group = group )

  if(length(tempNames) == 0){
    message("All names are correct")
    return(x)
  }
  else{
    wrong <- tempNames$wrongNames
    correct <- tempNames$correctNames
    result <- as.character(x$Taxa)
    n <- length(wrong)
    if ( length( wrong ) != length( correct ) ) {
      stop("pattern and replacement do not have the same length.")
    }
    for(i in 1:n){
      result[which(result == wrong[i])] <- correct[i]
    }

    x$Taxa <- result
    if( any( x$Taxa == "REMOVe" ) ){
      x <- x[ x$Taxa != "REMOVe" , ]
      }

    x$Taxa <- factor( sub( "_" , " ", x$Taxa ) )

    # write.table if needed by the user, allowing to overwrite or not
    if( write.table == TRUE ){
      dfNameChange <- paste0( dfName , "_mod" , ".txt" )
      if( dfNameChange %in% list.files(path = ".") ){
        choice.list <-  c( "Overwrite" , "Quit" )
        res <- menu( choice.list , title = paste( "Overwrite the existing " , dfNameChange , " file?" , sep = ""  ) )
        if( res == 1) { write.table( x , file = dfNameChange , quote = TRUE , append = FALSE , row.names = FALSE ) }
        else { stop( "Operation stopped by the user" ) }
      } else{ write.table( x , file = dfNameChange , quote = TRUE , append = FALSE , row.names = FALSE ) }
    }

    x
  }
}
