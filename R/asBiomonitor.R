#' asBiomonitor
#'
#' This function merge the user dataframe with the reference database. Options to improve or replace the reference database are provided.
#' @param x a data.frame as specified in details
#' @param group biotic group of interest. Possible values are "mi" for macroinvertebrates and "mf" for macrophytes. It overwrite = F group will not be considered
#' @param dref allow the user to improve (if overwrite = F) or to replace (if overwrite = T) the reference database
#' @param overwrite if set to T replace the reference database with the one provided by the user
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' asBiomonitor checks the correctness of taxa names in the data.frame provided by the user. If names are correct the function will process the data.frame to a biomonitor object, otherwise it provide suggestion for correct names. If dfref = T a custom dictionary will be saved in the working directory.
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{quickRename}}
#' @examples
#' data(macro_ex)
#' asBiomonitor(macro_ex, group = "mi")


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
    ref <- rbind(ref, dfref)
  }

  # allow the user to update the database replacing the reference database with is own reference database
  if(is.null(dfref) == F & overwrite == T){
    ref <- dfref
  }

  x <- aggregate(. ~ Taxa, x, FUN = sum)
  userTaxa <- x$Taxa

  # create a new dictionary to be used when user add taxa to the database


  # change the name of taxa to lowercase and capital letter
  userTaxaCap <- sapply(userTaxa, capWords, USE.NAMES = F)

  # changes various flavours of Hydracarina to Trombidiformes
  hydrac <- c("Hydracarina", "Hydracnidia", "Acariformes")
  hydrac_temp <- userTaxaCap %in% hydrac
  if(length(which(hydrac_temp == T)) != 0){
    userTaxaCap[which(hydrac_temp)] <- "Trombidiformes"
  }

  x$Taxa <- userTaxaCap
  if(is.null(dfref) == T){
    x <- rename(x, groups = group)
  }

  if(is.null(dfref) == F & overwrite == F){
    newDictio(ref)
    x <- rename(x, customx = T, groups = group)
  }

  if(is.null(dfref) == F & overwrite == T) {
    newDictio(ref)
    x <- rename(x, customx = T)
  }

  taxa_def <- merge(ref, x, by = "Taxa", all = F)

  if(group == "mi" & overwrite == F){
    class(taxa_def) <- c("biomonitoR", "mi")
  }
  if(group == "mf" & overwrite == F){
    class(taxa_def) <- c("biomonitoR", "mf")
  }
  if(overwrite == T){
    class(taxa_def) <- c("biomonitoR", "custom")
  }



  if(length(which(hydrac_temp == T)) != 0 ){
    message("Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes")
    taxa_def
  }

  else{ taxa_def }
}
