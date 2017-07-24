#' aggregatoR
#'
#' This function prepares data for further calculations.
#' @param x results of function asBiomonitoR
#' @keywords aggregatoR
#' @details
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)

aggregatoR <- function(x){
  if(class(x)!="biomonitoR"){
    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  tx <- c("Class", "Order", "Family", "Genus", "Species", "Taxa")
  stz <- x[!(names(x) %in% tx)]
  stz_n <- names(stz)     # station names

  # Counting Classes
  cla.agg <- aggregate(stz, by=list(x$Class),sum)
  levels(cla.agg$Group.1)[levels(cla.agg$Group.1)==""] <- "unassigned"
  names(cla.agg) <- c("Class",stz_n)

  # Counting Orders
  ord.agg <- aggregate(stz, by=list(x$Order),sum)
  levels(ord.agg$Group.1)[levels(ord.agg$Group.1)==""] <- "unassigned"
  names(ord.agg) <- c("Order",stz_n)

  # Counting Families
  fam.agg <- aggregate(stz, by=list(x$Family),sum)
  levels(fam.agg$Group.1)[levels(fam.agg$Group.1)==""] <- "unassigned"
  names(fam.agg) <- c("Family",stz_n)

  # Counting Genus
  gen.agg <- aggregate(stz, by=list(x$Genus),sum)
  levels(gen.agg$Group.1)[levels(gen.agg$Group.1)==""] <- "unassigned"
  names(gen.agg) <- c("Genus",stz_n)

  # Counting Species
  spe.agg <- aggregate(stz, by=list(x$Species),sum)
  levels(spe.agg$Group.1)[levels(spe.agg$Group.1)==""] <- "unassigned"
  names(spe.agg) <- c("Species",stz_n)

  # Taxa
  tax.agg <- aggregate(stz, by=list(x$Taxa),sum)
  levels(tax.agg$Group.1)[levels(tax.agg$Group.1)==""] <- "unassigned"
  names(tax.agg) <- c("Taxon",stz_n)

  tree.agg <- aggregate(stz, by=list(x$Class,x$Order, x$Family, x$Genus, x$Species, x$Taxa),sum)
  names(tree.agg) <- c(tx,stz_n)

  temp <- list(cla.agg,ord.agg,fam.agg,gen.agg,spe.agg,tax.agg,tree.agg)
  names(temp) <- c( "Class", "Order", "Family", "Genus", "Species", "Taxa", "Tree")
  class(temp) <- "biomonitoR"
  return(temp)
}
