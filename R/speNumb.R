#' speNumb
#'
#' Functions for calculating species, genus, family and order Richness
#' @param x results of function aggregatoR
#' @keywords speNumb, genNumb, famNumb, ordNumb
#' @details By now only species, genus and family richness calculation are reliable. This is because order assignment for order in the reference database is not completely covered. Unassigned taxon are exluded from the calculations.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' genNumb(data.agR)

speNumb <- function(x){
  spe <- x[["Species"]]
  if("unassigned" %in% spe[,1]){
    z <- which(spe$Species=="unassigned")
    spe <- spe[-z,] # remove unassigned row from the species count
  }
  nspe <- apply(spe[,-1], 2, FUN=function(x){length(x[x>0])})
  return(nspe)
}


genNumb <- function(x){
  gen <- x[["Genus"]]
  if("unassigned" %in% gen[,1]){
    z <- which(gen$Genus=="unassigned")
    gen <- gen[-z,] # remove unassigned row from the species count
  }
  ngen <- apply(gen[,-1], 2, FUN=function(x){length(x[x>0])})
  return(ngen)
}

famNumb <- function(x){
  fam <- x[["Family"]]
  if("unassigned" %in% fam[,1]){
    z <- which(fam$Family=="unassigned")
    fam <- fam[-z,] # remove unassigned row from the species count
  }
  nfam <- apply(fam[,-1], 2, FUN=function(x){length(x[x>0])})
  return(nfam)
}


ordNumb <- function(x){
  ord <- x[["Order"]]
  if("unassigned" %in% ord[,1]){
    z <- which(ord$Order=="unassigned")
    ord <- ord[-z,] # remove unassigned row from the species count
  }
  nord <- apply(ord[,-1], 2, FUN=function(x){length(x[x>0])})
  return(nord)
}


abu <- function(x){
  nord <- apply(x[["Order"]][,-1], 2, FUN=sum)
  return(nord)
}
