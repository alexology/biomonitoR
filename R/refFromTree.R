#' refFromTree
#'
#' This function transforms a taxonomic tree to a reference database suitable for biomonitoR.
#' @param x taxonomic tree
#' @keywords refFromTree
#' @export
#' @examples
#' data(Tree)
#' ref_custom <- refFromTree(Tree)



refFromTree <- function(x){
  
  n <- ncol(x)
  
  if(n == 1){
    stop("data.frame with 1 column are not allowed!")
  }
  
  taxa.col <- c("Phylum",	"Class",	"Subclass",	"Order",	"Family",	"Subfamily",
                "Tribus",	"Genus",	"Species",	"Subspecies",	"Taxa")
  
  empty.df <- as.data.frame( matrix(ncol = length(taxa.col), nrow = 0  ) )
  colnames(empty.df) <- taxa.col
  df <- empty.df
  
  x[] <- lapply(x, as.character)
  # colnames of the reference database provided by the user
  cref.name <- colnames(x)
  
  # transform to capital letter
  cref.name <- sapply(cref.name, capWords, USE.NAMES = F)
  
  # test if cref.name are different from colnames accepted by biomonitoR
  if(sum(! cref.name %in% taxa.col ) != 0 ){
    stop("Provide valid column names")
  }

  for(i in 1:n){
    temp.name <- colnames( x[ , i, drop = F] )
    temp.pos <- which(cref.name == temp.name) 
    temp <- x [ , 1:temp.pos, drop = F]
    # remove rows with empty cells
    temp <- temp[ which(temp[ , temp.name] != ""), , drop = F]
    temp.un <- unique(temp)
    empty.df <- merge(empty.df, temp.un, all.y  = T, sort = F)
    empty.df$Taxa <- temp.un[ , temp.name]
    df <- rbind.data.frame(df, empty.df)
  }
  
  df[is.na(df)] <- ""
  df <- as.data.frame( unclass( df ) )
  
  # remove leading and final spaces
  df <- sapply(df, trim, USE.NAMES = F)
  
  # remove NA originating from empty columns
  if(length(df[is.na(df)]) > 0){
    df[is.na(df)] <- ""
  }
  
  # check for duplicates or errors
  
  df <- as.data.frame( df )
  
  s.mes <- checkTree(df)
  if(is.null(s.mes) == T){
    stop(mes)
  }
  return( df )
}
