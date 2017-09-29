#' asBiomonitor
#'
#' This function merge the user dataframe with the reference database. Options to improve or replace the reference database are provided.
#' @param x a data.frame as specified in details
#' @param dref allow the user to improve (if overwrite = F) or to replace (if overwrite = T) the reference database
#' @param overwrite if set to T replace the reference database with the one provided by the user
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' asBiomonitor check the correctness of taxa names in the data.frame provided by the user. If names are correct the function will process the data.frame to a biomonitor object, otherwise it provide suggestion for correct names. If dfref = T a custom dictionary will be save in the working directory. 
#' @export
#' @seealso \code{\link{rename}}
#' @examples
#' data(macro_ex)
#' asBiomonitor(macro_ex)


asBiomonitor <- function (x, dfref = NULL, overwrite = F ) 
{
  
  # check if user database contains taxa of the reference database
  if(is.null(dfref) == F){
    temp.ref <- ref$Taxa
    temp.dfref <- dfref$Taxa
    both <- temp.ref[which(temp.dfref %in% temp.ref)]
    if(length(both) > 0){
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      return("user reference database contains taxa of the default reference database, please consider overwrite = T")
    }
    
  }
  
  # allow the user to update the database adding taxa to reference database
  if(is.null(dfref) == F & overwrite == F){
    ref <- rbind(ref, dfref)
  } 
  
  # allow the user to update the database replacing the reference database with is own reference database
  if(is.null(dfref) == F & overwrite == T){
    uref <- dfref
  }
  
  x <- aggregate(. ~ Taxa, x, FUN = sum)
  userTaxa <- x$Taxa
  
  # create a new dictionary to be used when user add taxa to the database
  
  
  # change the name of taxa to lowercase and capital letter
  userTaxaCap <- sapply(userTaxa, capWords, USE.NAMES = F)
  
  # changes various flavours of Hydracarina to Trombidiformes
  hydrac <- c("Hydracarina", "Hydracnidia", "Acariformes")
  hydrac_temp <- userTaxaCap %in% hydrac
  if(length(which(hydrac_temp == T)) != 0 ){
    userTaxaCap[which(hydrac_temp)] <- "Trombidiformes"
  }
  
  x$Taxa <- userTaxaCap
  if(is.null(dfref) == T){
    x <- rename(x)
  } else {
    if(is.null(dfref) == F & overwrite == F){
      newDictio(ref)
      x <- rename(x, custom = T)
    } else {
      newDictio(uref)
      x <- rename(x, custom = T)
    }
    
  }
  
  temp <- merge(ref, x, by = "Taxa", all = F)
  temp_valid <- temp[which(temp$Taxonomic_Status=="yes"),]
  temp_novalid <- temp[which(temp$Taxonomic_Status=="no"),]
  taxa_def <- as.list(temp_valid[,-which(names(temp) %in% c("Taxonomic_Status"))])
  tx <- c("Class",  "Subclass", "Order", "Family", "Genus", "Species", "Taxonomic_Status")
  taxa_def$novalid <- temp_novalid[!(names(temp_novalid) %in% tx)]
  
  class(taxa_def) <- "biomonitoR"
  
  if(length(which(hydrac_temp == T)) != 0 ){
    message("Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes")
    taxa_def
  }
  
  else{ taxa_def }
}
