#' igold
#'
#' This function calculates the 1 - GOLD metric, where GOLD stands for Gastropoda, Oligochaeta and Diptera. This metric should decrease with increasing organic pollution (Pinto et al., 2004).
#' @param x results of function aggregatoR
#' @keywords ept
#' @details The metric 1 - GOLD is calculated as 1 minus the relative abundance of Gastropoda, Oligochaeta and Diptera. If a custom database is provided (see \code{\link{aggregatoR}}) please be sure that Gastropoda and Oligochaeta are submitted as Class and Diptera as Order, otherwise the gold calculation will be meaningless.
#' @importFrom stats aggregate
#' @references Pinto, P., Rosado, J., Morais, M., & Antunes, I. (2004). Assessment methodology for southern siliceous basins in Portugal. In Integrated Assessment of Running Waters in Europe (pp. 191-214). Springer, Dordrecht.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' igold(data.agR)


igold <- function (x){

  # check if the object x is of class "biomonitoR"
  classCheck(x, group = "mi")


  x_gold <- x[["Tree"]]
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz <- x_gold[!(names(x_gold) %in% tx)]
  stz_n <- names(stz)     # station names
  gold_taxa <- x_gold[which(x_gold$Class == "Gastropoda" |
                            x_gold$Subclass == "Oligochaeta" | x_gold$Order == "Diptera"), , drop=F]
  if(nrow(gold_taxa)==0){
    gold_taxa[1,-1] <- rep(0, ncol(gold_taxa)-1)
    return(gold_taxa[,-1])
  } else {
    gold_temp <- gold_taxa[,c("Taxa" , stz_n)]
    colnames(gold_temp)[1] <- "selection"
    levels(gold_temp$selection)[levels(gold_temp$selection)==""] <- "unassigned"
    gold_temp.agg <- aggregate(. ~ selection, gold_temp, FUN=sum)
    if("unassigned" %in% gold_temp.agg[,1]){
      z <- which(gold_temp.agg[,1]=="unassigned")
      gold_temp.agg<- gold_temp.agg[-z,] # remove unassigned row from the species count
    }
    temp <- 1 - apply(gold_temp.agg[ , -1 , drop = FALSE ], 2 ,sum) / abu(x)
    return( temp )
  }
}
