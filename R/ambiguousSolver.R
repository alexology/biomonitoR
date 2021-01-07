#' ambiguousSolver
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("deprecated") }
#'
#' Solve the ambiguous assignment of taxa.
#'
#' @param x result of the function aggregatoR.
#' @param method only methods `RPKC` and `MCWP` are currently implemented in `biomonitoR`. See details.
#' @details Taxonomic dataset often contains ambiguous taxa due to the difficulites in identify damaged and juvenile organisms.
#' For instance Chironominae and Chironomidae could be reported in the same sample. Several techinques has been
#' proposed to solve this issue (see Cuffney et al., 2007).
#' Currently ambiguousSolver two only one algorithms to solve the problem of ambiguous assignments.
#' `RPKC` means Remove Parents - Keep Children. It works by removing all the higher taxonomic levels while keeping children.
#' `MCWP` means Merge Children with Parents. It works by assignign the abundances of lower taxonomic levels to those of higher taxonomic levels.
#' @keywords aggregate_taxa
#' @references Cuffney, T. F., Bilger, M. D., & Haigler, A. M. (2007). Ambiguous taxa: effects on the characterization and interpretation of invertebrate assemblages. Journal of the North American Benthological Society, 26(2), 286-307.
#' @export
#' @seealso \code{\link{aggregatoR}}


ambiguousSolver <- function( x , method = "MCWP" ){

  .Deprecated("solve_ambiguous")

  # check if the object d is of class "biomonitoR"
  classCheck( x )

  # Extract the taxonomic tree
  tree <- x[[ "Tree" ]]

  # order tree from the rows containing species to the rows containing only Phylum ( if present, obviously )
  # order also alphabetically ( maybe not needed )
  n.tl <- apply( tree[ , 1:10 ] , 1 , function( x ) max( which( x != "" ) ) )
  tree <- tree[ order( n.tl , decreasing = TRUE ) , ]

  # Extract the taxa list as charcater and get the position of the column called taxa in the taxonomic tree
  taxa <- as.character( tree$Taxa )
  taxa.pos <- which( names( tree ) == "Taxa" )

  # intialize the vector that will store the taxa already processed
  # this should speed up the computation time
  taxa.acc <- c()

  # RPKC: basically, for each taxon in the taxonomic list we extract the taxonomic tree from taxlev -1 to Phylum
  # taxlev -1 represents the taxonomic level before the taxonomic level of the working taxon: temp <- temp[ -length( temp ) ]
  # if any of the taxa along the tree from taxlev -1 until Phylum is present in the column taxa we have a problem of
  # ambiguous taxa and we store the results in the taxa.acc vector. The higher taxonomic levels are then deleted from the database.

  if( method == "RPKC" ){

    for( i in 1:length( taxa ) ){
      temp <- tree[  i , 1:(taxa.pos - 1) ]
      temp <- temp[ temp != "" ]
      temp <- temp[ -length( temp ) ]
      if( any( temp %in% taxa ) ){
        taxa.acc <- c( taxa.acc , temp[ temp %in% taxa ] )
      }
    }

    # get unique higher taxa names
    taxa.acc <- unique( taxa.acc )

    # make a vector of the same length as the column taxa, set value to TRUE
    # and set ambiguous taxa to FALSE. The orignal taxa list is used
    # instead of the ordered
    taxa.pc <- rep( TRUE, length( taxa ) )
    taxa.pc[ x[[ "Taxa" ]][ , "Taxon" ] %in% taxa.acc ] <- FALSE

    # parents are thus removed from the database
    DF <- x[[ "Taxa" ]][ taxa.pc , ]
    names( DF )[ 1 ] <- "Taxa"
  }

  # MCWP: basically, the taxonomic tree of each taxon is retrieved. The presence of ambiguos taxa in the taxonomic tree is
  # checked against the Taxa column of the user database. If present we have ambiguous taxa to solve, otherwise no.
  # The taxon name of the i-th row is set to the highest taxonomic level if ambiguous parents are present, otherwise no.


  if( method == "MCWP" ){
    for( i in 1:length( taxa ) ){
      if( any( taxa[ i ] %in% taxa.acc ) ){
        next
      } else {
        temp <- tree[ tree$Taxa %in% taxa[ i ] , -c( taxa.pos : ncol( tree ) ) ]
        temp.tax_lev <- names( temp[ , temp != "" ] )
        temp <- temp[ temp != "" ]
        temp <- temp[ -length( temp ) ]
        if( ! any( temp %in% taxa ) ){
          taxa.acc <- c( taxa.acc , taxa[ i ] )
          next
        } else {
          t.taxa <- min( which( temp %in% taxa ) )
          taxa.up <- temp[ t.taxa ]
          tax_lev.up <- temp.tax_lev[ t.taxa ]
          taxa[ tree[ , tax_lev.up ] %in% taxa.up ] <- taxa.up
          taxa.acc <- c( taxa.acc , taxa.up )
        }
      }
    }

    # build the data.frame and order it as the original
    DF <-  data.frame( Taxa = as.factor( taxa ) , tree[ , ( taxa.pos + 1 ):ncol( tree ) , drop = FALSE ] )
    DF <- DF[ match( 1:nrow( tree ) , as.numeric( rownames( DF ) ) ) , ]

  }

  # aggregate to avoid duplicated rownames
  DF <- aggregate( . ~ Taxa , data = DF , FUN = sum )

  if( inherits( x , "bin" ) ){
    DF <- to_bin( DF )
  }

  DF
}
