#' asBiomonitor
#'
#' This function prepares data for further calculations.
#' @param x a data.frame as specified in details
#' @keywords asBiomonitor
#' @details data.frame must have a column called "Taxa" where put species, genus or family names. See data(macro_ex) for an example dataset.\cr
#' asBiomonitor check the correctness of taxa names in the data.frame provided by the user. If names are correct the function will process the data.frame to a biomonitor object, else it suggest to rename the wrong names with the function \code{\link{rename}}.
#' @export
#' @seealso \code{\link{rename}}
#' @examples
#' data(macro_ex)
#' asBiomonitor(macro_ex)


asBiomonitor <- function (x, dfref = NULL, overwrite = F ) 
{
  
  # allow the user to update the database adding taxa to reference condition
  if(is.null(dfref) == F & overwrite == F){
    ref <- rbind(ref, dfref)
  } 
  
  # allow the user to update the database replacing the reference database with is own reference database
  if(is.null(dfref) == F & overwrite == T){
    ref <- dfref
  }
  
  x <- aggregate(. ~ Taxa, x, FUN = sum)
  userTaxa <- x$Taxa
  
  # cahnge the name of taxa to lowercase and capital letter
  userTaxaCap <- sapply(userTaxa, capWords, USE.NAMES = F)
  
  # changes various flavours of Hydracarina to Trombidiformes
  hydrac <- c("Hydracarina", "Hydracnidia", "Acariformes")
  hydrac_temp <- userTaxaCap %in% hydrac
  if(length(which(hydrac_temp == T)) != 0 ){
    userTaxaCap[which(hydrac_temp)] <- "Trombidiformes"
  }
    
  x$Taxa <- userTaxaCap
  x <- rename(x)
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
