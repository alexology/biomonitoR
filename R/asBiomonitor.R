#' asBiomonitor
#'
#' This function merge the user dataframe with the reference database (be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015)). Options to improve or replace the reference database are provided.
#'
#' @param x a data.frame as specified in details
#' @param group biotic group of interest. Possible values are "mi" for macroinvertebrates and "mf" for macrophytes. This option will not be considered if overwrite i set to FALSE.
#' @param dfref allow the user to improve (if overwrite = F) or to replace (if overwrite = T) the reference database.
#' @param overwrite if set to T the reference database is replaced with the one provided by the user.
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' asBiomonitor checks the correctness of taxa names in the data.frame provided by the user. If names are correct the function will process the data.frame to a biomonitor object, otherwise it will provide suggestions for correct names. If dfref = T a custom dictionary will be saved in the working directory.
#' When using the default database the taxa names provided by the user need to be consistent with the taxonomy of [freshwaterecology.info](https://www.freshwaterecology.info/), otherwise the user is asked to exit. This behaviour is to assure consistency with other functions implemented in biomonitoR.
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{quickRename}}
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www. freshwaterecology.
#'   info-An online tool that unifies, standardises and codifies more than
#'   20,000 European freshwater organisms and their ecological preferences.
#'   Ecological indicators, 53, 271-282.
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex, group = "mi")


asBiomonitor <- function (x, group = "mi", dfref = NULL, overwrite = F )
{
  # check if user database contains a column called Taxa
  if(!"Taxa" %in% names(x)){
    stop("Column called Taxa needed")
  }


  # check if columns other than Taxa are numeric
  # position of column Taxa
  col.taxa <- which(names(x) == "Taxa")
  col.class <- sapply(x[, -col.taxa ], is.numeric)
  if(any(col.class == FALSE)){
    stop("Non-numeric columns are not allowed")
  }

  # check if x is a presence absence data.frame
  check.pa <- any( x[ , -which( "Taxa" %in% colnames( x ) ) , drop = FALSE ]  > 1 )

  if(group == "mi"){
    ref <- mi_ref
  }

  if(group == "mf"){
    ref <- mf_ref
  }

  # check if user database contains taxa of the reference database
  if(is.null(dfref) == F & overwrite == F){
    temp.ref <- ref$Taxa
    temp.dfref <- dfref$Taxa
    both <- temp.ref[which(temp.dfref %in% temp.ref)]
    if(length(both) > 0){
      stop("user reference database contains taxa of the default reference database, please consider overwrite = T")
    }

  }

  # allow the user to update the database adding taxa to reference database
  if(is.null(dfref) == F & overwrite == F){
    dfref[ is.na( dfref ) ] <- ""
    dfref <- as.data.frame( unclass( dfref ) )
    ref <- rbind(ref, dfref)
  }

  # allow the user to update the database replacing the reference database with is own reference database
  if(is.null(dfref) == F & overwrite == T){
    dfref[ is.na( dfref ) ] <- ""
    dfref <- as.data.frame( unclass( dfref ) )
    ref <- dfref
  }

  x <- aggregate(. ~ Taxa, x, FUN = sum)
  userTaxa <- trimws(x$Taxa)

  # create a new dictionary to be used when user add taxa to the database


  # change the name of taxa to lowercase and capital letter
  userTaxaCap <- sapply(userTaxa, capWords, USE.NAMES = F)

  # initialize message for Trombidiformes
  mes <- NULL

  if(group == "mi"){
    # changes various flavours of Hydracarina to Trombidiformes
    hydrac <- c("Hydracarina", "Hydracnidia", "Acariformes")
    hydrac_temp <- userTaxaCap %in% hydrac
    if(length(which(hydrac_temp == T)) != 0){
      userTaxaCap[which(hydrac_temp)] <- "Trombidiformes"
      mes <- "Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes"
    }
  }


  x$Taxa <- userTaxaCap
  if(is.null(dfref) == T){
    x <- rename(x, group = group)
  }

  if(is.null(dfref) == F & overwrite == F){
    newDictio(ref)
    x <- rename(x, customx = T, group = group)
  }

  if(is.null(dfref) == F & overwrite == T) {
    newDictio(ref)
    x <- rename(x, customx = T)
  }

  # aggregate another time to take into account name changes
  x <- aggregate(. ~ Taxa, x, FUN = sum)
  if( ! check.pa ) ( x <- data.frame( x[ , 1 , drop = FALSE], ( x[ , -1, drop = FALSE ] > 0 ) * 1 ) )

  taxa_def <- merge(ref, x, by = "Taxa", all = F)

  if(group == "mi" & is.null( dfref ) == TRUE){
    class(taxa_def) <- c("biomonitoR", "mi")
  }
  if(group == "mf" & is.null( dfref ) == TRUE){
    class(taxa_def) <- c("biomonitoR", "mf")
  }
  if(is.null( dfref ) == FALSE){
    class(taxa_def) <- c("biomonitoR", "custom")
  }

  if( is.null(mes) == F ){
    message( mes )
    taxa_def
  }

  else{ taxa_def }
}
