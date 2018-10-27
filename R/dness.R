#' dness
#'
#' This function calculates the taxonomic diversity, taxonomic distinctness and variation of taxonomic distinctness.
#' @param x results of function aggregatoR
#' @param complete if TRUE unassigned taxon are removed from the taxonomic list of the desired taxonomic level set with taxLev
#' @param taxLev taxonomic level from which the indices are to be calculated
#' @param method available methods are delta (taxonomic diversity), delta.st (taxonomic distinctness) and delta.bin (variation in taxonomic distinctness)
#' @param taxatree if TRUE dness return a list with results and taxonomic tree
#' @keywords ept
#' @references Clarke, K. R., & Warwick, R. M. (1998). A taxonomic distinctness index and its statistical properties. Journal of applied ecology, 35(4), 523-531.
#' @references Clarke, K. R., & Warwick, R. M. (2001). A further biodiversity index applicable to species lists: variation in taxonomic distinctness. Marine ecology Progress series, 216, 265-278.
#' @details see Clarke and Warwick (1998) an Clarke and Warwick (2001).
#' @note USE WITH CAUTION, STILL IN DEVELOPMENT.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' dness(data.agR, complete = TRUE)




dness <- function(x, complete = FALSE, taxLev = "Genus", method = "delta", taxatree = FALSE) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # check if there are unassigned taxon in the taxonomic tree
  # if present stop the function, unless the user specify unassigned = T

  df <- x[[taxLev]]
  
  # check if the selected taxonomic level contains a least one taxa
  if( nrow( df ) == 1 & df[ 1 , 1 ] == "unassigned" ){ stop("The selected taxonomic level is empty") }
  if( taxLev == "Taxa"){ stop("taxLev cannot be equal to Taxa") }
  
  
  st.names <- names(df[, -which(names(df) %in% c(taxLev)), drop = FALSE])   # sites name
  tax <- df[ , taxLev]    # taxa names at the desired taxonomic level taxLev
  che.ck <- "unassigned" %in% tax   # check if unassigned is present at the desired taxonomic level
  tree.c <- x[["Tree"]]
  # exclude taxonomic levels with missing information. This feature could be improved in the next future
  exclude <- -which(names(tree.c) %in% c("Taxa", "Subclass", "Subfamily", "Tribus", "Subspecies"))
  tree <- tree.c[ , exclude, drop = FALSE]
  # taxonomic levels to retain
  tax.pos <- names(tree[, 1:which(names(tree) == taxLev)])
  # ref.c <- tree.c[ , c(which(names(tree.c) %in%  c( tax.pos, "Taxa" ) ))]
  ref.c <- tree.c[ , c(which(names(tree.c) %in%  tax.pos ))]
  # remove duplicated rows
  ref.c <- ref.c[ ! duplicated( ref.c ) ,  ]
  
  if(complete == FALSE & che.ck == TRUE) {stop("Unassigned taxon at the desired taxonomic level")}
  if(complete == FALSE & che.ck == FALSE){
    # retain only the taxonomic levels from order to taxLev and sites
    taxtree <- merge(ref.c, df, by = taxLev, all = FALSE)
    taxtree <- taxtree[, c( tax.pos , st.names ), drop = FALSE]
  }

  if(complete == TRUE & che.ck == TRUE){
    df <- df[ -which( df[ , taxLev] == "unassigned" ) , ]
    taxtree <- merge(ref.c, df, by = taxLev, all = FALSE)
    taxtree <- taxtree[, c( tax.pos , st.names ), drop = FALSE]
  }
  
  if(complete == TRUE & che.ck == FALSE){
    taxtree <- merge(ref.c, df, by = taxLev, all = FALSE)
    taxtree <- taxtree[, c( tax.pos , st.names ), drop = FALSE]
  }

  # remove empty columns (maybe is not a necessary step)
  temp <- taxtree[ , which(names(taxtree) %in% c(tax.pos)), drop = FALSE]
  temptf <- temp != ""
  # remove columns with at least 1 empty cell
  check.col <- apply(temptf, 2, sum) / nrow(temptf)
  if(sum(check.col) != length(temp)){
    taxtree <- taxtree[, - which(check.col != 1), drop = FALSE]
  }
  
  # remove columns distinct for all species or with the same taxa name for all the rows (excluding taxLev)
  # -1 is to eclude taxLev from this computation
  tax.pos2 <- 1:(which(names(taxtree) == taxLev) - 1)
  tax.uniq <- apply( taxtree[ , tax.pos2 , drop = FALSE] , 2 , function(x) ( length(unique(x) ) ) )
  col.rm <- -which(tax.uniq == 1 | tax.uniq == nrow(taxtree) )
  if(length(col.rm) != 0) (taxtree <- taxtree[ , col.rm , drop = FALSE ])
  
  
  if(ncol(taxtree[, -which(names(taxtree) %in% c(st.names)), drop = FALSE]) == 1){
    stop("Reference database has not enough taxonomic levels to perform the analysis")
  }


  # calculating taxonomic distance
  tax.dis <- ddis(taxtree, st.names = st.names)
  sites <- taxtree[ , which(names(taxtree) %in% st.names), drop = FALSE]
  sites.bin <- sites
  sites.bin[sites.bin > 0 ] <- 1

  if(method == "delta") {res <- apply(sites, 2, function( x ) ( delta(x, dis = tax.dis)))}
  if(method == "delta.bin") {res <- apply(sites.bin, 2, function( x ) ( delta(x, dis = tax.dis)))}
  if(method == "delta.st") {res <- apply(sites, 2, function( x ) ( delta.st(x, dis = tax.dis)))}

  if(taxatree == FALSE) return(res)
  if(taxatree == TRUE) return(list(res, taxtree))
}
