#' ept
#'
#' This function calculates EPT richnes
#' @param x results of function aggregatoR
#' @keywords ept
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' ept(data.bio)


epti <- function (x, taxLev = "Family"){
  x_ept <- x[["Tree"]]
  tx <- c("Class", "Order", "Family", "Genus", "Species", "Taxa")
  stz <- x_ept[!(names(x_ept) %in% tx)]
  stz_n <- names(stz)     # station names
  ept_taxa <- x_ept[which(x_ept$Order == "Plecoptera" |
                            x_ept$Order == "Ephemeroptera" | x_ept$Order == "Trichoptera"),]
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
    temp <- colSums(ept_temp.agg[,-1,drop=F] > 1)
    return(temp)
  }
}
