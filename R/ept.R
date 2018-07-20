#' ept
#'
#' This function calculates the number of Ephemeroptera, Plecotera and Trichoptera (EPT) taxa at different taxonomic levels.
#' @param x results of function aggregatoR
#' @param taxLev the taxonomic level for calculating EPT richness.
#' @keywords ept
#' @details The parameter taxLev can be "Species", "Genus", "Family" or "Order".
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' ept(data.agR)


ept <- function (x, taxLev = "Family"){

  # check if the object d is of class "biomonitoR"

  # check if the object x is of class "biomonitoR"
  classCheck(x, group = "mi")


  x_ept <- x[["Tree"]]
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz <- x_ept[!(names(x_ept) %in% tx)]
  stz_n <- names(stz)     # station names
  ept_taxa <- x_ept[which(x_ept$Order == "Plecoptera" |
                            x_ept$Order == "Ephemeroptera" | x_ept$Order == "Trichoptera"), , drop=F]
  if(nrow(ept_taxa)==0){
    ept_taxa[1,-1] <- rep(0, ncol(ept_taxa)-1)
    return(ept_taxa[,-1])
  } else {
    ept_temp <- ept_taxa[,c(taxLev,stz_n)]
    colnames(ept_temp)[1] <- "selection"
    levels(ept_temp$selection)[levels(ept_temp$selection)==""] <- "unassigned"
    ept_temp.agg <- aggregate(. ~ selection, ept_temp, FUN=sum)
    if("unassigned" %in% ept_temp.agg[,1]){
      z <- which(ept_temp.agg[,1]=="unassigned")
      ept_temp.agg<- ept_temp.agg[-z,] # remove unassigned row from the species count
    }
    temp <- colSums(ept_temp.agg[,-1,drop=F] > 0)
    return(temp)
  }
}
