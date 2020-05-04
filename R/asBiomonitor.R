#' asBiomonitor
#'
#' This function merge the user dataframe with a reference database and suggest corrections for mispelled names.
#'
#' @param x a data.frame with a column called "Taxa" where store taxa names and samples on the other columns (see the example macro_ex).
#' @param group biotic group of interest. Possible values are `mi` for macroinvertebrates and `mf` for macrophytes. The choice will set the right reference database for the specified group.
#' This option will not be considered if a custom reference database is provided. Default to `mi`.
#' @param dfref a custom reference database that replaces the reference database.
#' @param FUN the function to be applied for aggregating rows with duplicated taxa names.Default to `sum`.
#' @keywords asBiomonitor
#' @details The function `asBiomonitor` checks the taxonomy of the data.frame provided by the user and suggests correction for mispelled names.
#' If one or more taxa names of `x` are not present in the reference database or the spell checker is not able to find any suggestion the user is asked to exit.
#' This behaviour is to assure consistency with other functions implemented in biomonitoR.
#' Default references databases are provided for macroinvertebrates and macrophytes.
#' Both databases heavily rely on the information provided by the [freshwaterecology.info](https://www.freshwaterecology.info/) website.
#' If `dfref` is not NULL a custom dictionary will be saved in the working directory to let the `asBiomonitor` function work correctly.
#' If you are unable to build a reference database by your own please check the function \code{\link{refFromTree}} for a possible solution.
#' `asBiomonitor`  returns an object of class `biomonitoR` togheter with one of the classes `mi`, `mf` or `custom` depending from the `group`
#' and `dfref` specifications. If the reference database is used, the class is `mi` for macroinvertebrates and `mf` for macrophytes.
#' If a reference database is provided the class is set to `custom`. The function \code{\link{quickRename}} works as the `asBiomonitor` but returns
#' a data.frame without the biomonitoR format.
#' `asBiomonitor` aggregates all the rows with the same name with the option `FUN` and converts all the `NA` to 0.
#' If only 1 and 0 are present `x` will be imported as presence-absence.
#' When `group = mi` Hydracarina, Hydracnidia or Acariformes are changed to Trombidiformes given the uncertain taxonomic status of this group.
#'
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{quickRename}} \code{\link{refFromTree}} \code{\link{quickRename}}
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www.freshwaterecology.info -
#' An online tool that unifies, standardises and codifies more than
#' 20,000 European freshwater organisms and their ecological preferences.
#' Ecological indicators, 53, 271-282.
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex, group = "mi")


asBiomonitor <- function ( x , group = "mi" , dfref = NULL , FUN = sum ){

  # check if user database contains a column called Taxa
  if( ! "Taxa" %in% names( x ) ){
    stop( "A column called Taxa is needed")
  }

  # check if columns other than Taxa are numeric
  # position of column Taxa
  col.taxa <- which( names( x ) == "Taxa" )
  col.class <- sapply( x[ , -col.taxa ], is.numeric )
  if( any( col.class == FALSE ) ){
    stop( "Non-numeric columns other than Taxa are not allowed" )
  }

  # check if x is a presence absence data.frame
  check.pa <- any( x[ , -which( "Taxa" %in% colnames( x ) ) , drop = FALSE ]  != 1 )

  # check if NAs are present. If present they are changed to 0.
  check.na <- any(  is.na( x ) )
  x[ is.na( x ) ] <- 0

  # set the reference database for the specified group
  if( group == "mi" ){
    ref <- mi_ref
  }

  if(group == "mf"){
    ref <- mf_ref
  }

  # allow the users to use their own reference database
  if( ! is.null( dfref ) ){
    dfref[ is.na( dfref ) ] <- ""
    dfref <- as.data.frame( unclass( dfref ) )
    ref <- dfref
    group <- "custom"
  }

  # remove  leading and/or trailing whitespaces
  userTaxa <- trimws( x$Taxa )

  # change the name of taxa to lowercase and capital letter
  userTaxaCap <- sapply( userTaxa , capWords , USE.NAMES = F )

  # initialize message for Trombidiformes for the two checks
  mes <- NULL
  mes2 <- NULL

  # change the Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes
  if( group == "mi" ){
    # changes various flavours of Hydracarina to Trombidiformes
    hydrac <- c( "Hydracarina" , "Hydracnidia" , "Acariformes" )
    hydrac_temp <- userTaxaCap %in% hydrac
    if( length( which( hydrac_temp == T ) ) != 0 ){
      userTaxaCap[ which( hydrac_temp ) ] <- "Trombidiformes"
      mes <- "Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes"
    }
  }

  x$Taxa <- userTaxaCap

  # if dfref is NULL use the rename function to check for mispelled names and suggest for correct names
  if( is.null( dfref ) ){
    x <- rename( x , group = group )
  }

  # if dfref is not NULL create the dictionary in the user working directory
  # and use the rename function to check for mispelled names and suggest for correct names
  if( ! is.null( dfref ) ) {
    newDictio( ref )
    x <- rename( x , customx = T )
  }

  # check for Hydracarina again, in case the user mispelled the name
  userTaxaCap <- as.character( x$Taxa )

  if( group == "mi"  ){
    # changes various flavours of Hydracarina to Trombidiformes
    hydrac <- c( "Hydracarina" , "Hydracnidia" , "Acariformes" )
    hydrac_temp <- userTaxaCap %in% hydrac
    if( length( which( hydrac_temp == T ) ) != 0 ){
      userTaxaCap[ which( hydrac_temp ) ] <- "Trombidiformes"
      mes2 <- "Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes"
    }
  }

  x$Taxa <- userTaxaCap

  # aggregate once again to take into account name changes
  x <- aggregate( . ~ Taxa, x , FUN = FUN )
  if( ! check.pa ){
    x <- data.frame( x[ , 1 , drop = FALSE ], ( x[ , -1, drop = FALSE ] > 0 ) * 1 )
    message( "data imported as presence absence" )
  }

  if( isTRUE( check.na ) ){
    message( "NA detected, transformed to 0" )
  }

  # merge reference database to the user data.frame
  taxa_def <- merge( ref, x , by = "Taxa" , all = F )

  # assign classes to the data.frame
  if(group == "mi" & is.null( dfref ) == TRUE){
    class(taxa_def) <- c("biomonitoR", "mi" )
  }
  if(group == "mf" & is.null( dfref ) == TRUE){
    class(taxa_def) <- c( "biomonitoR" , "mf" )
  }
  if(is.null( dfref ) == FALSE){
    class(taxa_def) <- c( "biomonitoR", "custom" )
  }

  if( ( ! is.null( mes ) ) | ( ! is.null( mes2 ) )  ){
    message( "Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes" )
    taxa_def
  } else{ taxa_def }

}
