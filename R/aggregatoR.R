#' aggregatoR
#'
#' This function prepares data for further calculations.
#' @param z results of function asBiomonitoR
#' @keywords aggregatoR
#' @details
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)

aggregatoR <- function (z)
{
  if (class(z) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  x <- z
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz <- x[!(names(x) %in% tx)]
  stz_n <- names(stz)
  phy.agg <- aggregate(stz, by = list(x$Phylum), sum)
  levels(phy.agg$Group.1)[levels(phy.agg$Group.1) == ""] <- "unassigned"
  names(phy.agg) <- c("Phylum", stz_n)
  cla.agg <- aggregate(stz, by = list(x$Class), sum)
  levels(cla.agg$Group.1)[levels(cla.agg$Group.1) == ""] <- "unassigned"
  names(cla.agg) <- c("Class", stz_n)
  scla.agg <- aggregate(stz, by = list(x$Subclass), sum)
  levels(scla.agg$Group.1)[levels(scla.agg$Group.1) == ""] <- "unassigned"
  names(scla.agg) <- c("Subclass", stz_n)
  ord.agg <- aggregate(stz, by = list(x$Order), sum)
  levels(ord.agg$Group.1)[levels(ord.agg$Group.1) == ""] <- "unassigned"
  names(ord.agg) <- c("Order", stz_n)
  fam.agg <- aggregate(stz, by = list(x$Family), sum)
  levels(fam.agg$Group.1)[levels(fam.agg$Group.1) == ""] <- "unassigned"
  names(fam.agg) <- c("Family", stz_n)
  sfam.agg <- aggregate(stz, by = list(x$Subfamily), sum)
  levels(sfam.agg$Group.1)[levels(sfam.agg$Group.1) == ""] <- "unassigned"
  names(sfam.agg) <- c("Subfamily", stz_n)
  tri.agg <- aggregate(stz, by = list(x$Tribus), sum)
  levels(tri.agg$Group.1)[levels(tri.agg$Group.1) == ""] <- "unassigned"
  names(tri.agg) <- c("Tribus", stz_n)  
  gen.agg <- aggregate(stz, by = list(x$Genus), sum)
  levels(gen.agg$Group.1)[levels(gen.agg$Group.1) == ""] <- "unassigned"
  names(gen.agg) <- c("Genus", stz_n)
  spe.agg <- aggregate(stz, by = list(x$Species), sum)
  levels(spe.agg$Group.1)[levels(spe.agg$Group.1) == ""] <- "unassigned"
  names(spe.agg) <- c("Species", stz_n)
  sspe.agg <- aggregate(stz, by = list(x$Subspecies), sum)
  levels(sspe.agg$Group.1)[levels(sspe.agg$Group.1) == ""] <- "unassigned"
  names(sspe.agg) <- c("Subspecies", stz_n)
  tax.agg <- aggregate(stz, by = list(x$Taxa), sum)
  levels(tax.agg$Group.1)[levels(tax.agg$Group.1) == ""] <- "unassigned"
  names(tax.agg) <- c("Taxon", stz_n)
  tree.agg <- aggregate(stz, by = list(x$Phylum, x$Class, x$Subclass, x$Order, x$Family, x$Subfamily,
                                       x$Tribus, x$Genus, x$Species, x$Subspecies, x$Taxa), sum)
  names(tree.agg) <- c(tx, stz_n)
  temp <- list(phy.agg, cla.agg, scla.agg, ord.agg, fam.agg, sfam.agg, tri.agg, gen.agg, spe.agg,
               sspe.agg, tax.agg, tree.agg)
  names(temp) <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa", "Tree")
  class(temp) <- "biomonitoR"
  return(temp)
}
