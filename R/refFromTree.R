#' refFromTree
#'
#' This function transforms a taxonomic tree to a reference database suitable for biomonitoR.
#' @param x taxonomic tree. See \code{\link{Tree}} for an example.
#' @param group merge the user database with biomonitoR reference databases, default to "none". If duplicated Taxa names are present this function keeps the name provided by the user. Check the reliability of results when usin group = "mi" or group = "mf"
#' @keywords refFromTree
#' @export
#' @examples
#' data(Tree)
#' ref_custom <- refFromTree(Tree)



refFromTree <- function(x, group = "none"){

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
  cref.name <- sapply(cref.name, capWords, USE.NAMES = FALSE)

  # test if cref.name are different from colnames accepted by biomonitoR
  if(sum(! cref.name %in% taxa.col ) != 0 ){
    stop("Provide valid column names")
  }

  # reorder user data.frame according to biomonitoR taxa tree
  value.match <- match( taxa.col , names( x ) )
  x <- x[ , value.match[ ! is.na( value.match ) ] ]

  cref.name <- colnames( x )

  for(i in 1:n){
    temp.name <- colnames( x[ , i, drop = FALSE] )
    temp.pos <- which(cref.name == temp.name)
    temp <- x [ , 1:temp.pos, drop = FALSE]
    # remove rows with empty cells
    temp <- temp[ which(temp[ , temp.name] != ""), , drop = FALSE]
    temp.un <- unique(temp)
    empty.df <- merge(empty.df, temp.un, all.y  = TRUE, sort = FALSE)
    empty.df$Taxa <- temp.un[ , temp.name]
    df <- rbind.data.frame(df, empty.df)
  }

  df[is.na(df)] <- ""
  df <- as.data.frame( unclass( df ) )

  # remove leading and final spaces
  df <- sapply(df, trim, USE.NAMES = FALSE)

  # remove NA originating from empty columns
  if(length(df[is.na(df)]) > 0){
    df[is.na(df)] <- ""
  }

  # check for duplicates or errors

  df <- as.data.frame( df )
  df <- df[ ! duplicated( df ) , ]

  s.mes <- checkTree(df)
  if(is.null(s.mes) == FALSE){
    stop(s.mes)
  }

  if(group == "mi"){
    df <- rbind(df, mi_ref)
    df <- df[ !duplicated(df$Taxa) ,]
  }

  if(group == "mf"){
    df <- rbind(df, mf_ref)
    df <- df[ !duplicated(df$Taxa) ,]
  }

  return( df )

}
