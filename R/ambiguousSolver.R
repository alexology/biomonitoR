#' ambiguousSolver
#'
#' Solve the ambiguous assignment of taxa.
#'
#' @param x result of function aggregatoR.
#' @param method only method *RPKC* and *MCWP* is currently implemented. See details.
#' @details Taxonomic dataset often contains ambiguous taxa due to the difficulites in identify damaged and juvenile organisms. For instance Chironominae and Chironomidae could reported as in the same sample. Several techinques has been
#' proposed to solve this issue (see Cuffney et al., 2007). Currently ambiguousSolver has only one algorithm that allow to assing the abundances of lowr taxonomic level to those of higher taxonomic level.
#' @keywords aggregatoR
#' @references Cuffney, T. F., Bilger, M. D., & Haigler, A. M. (2007). Ambiguous taxa: effects on the characterization and interpretation of invertebrate assemblages. Journal of the North American Benthological Society, 26(2), 286-307.
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data( macro_ex )
#' data.bio <- asBiomonitor( macro_ex )
#' data.agR <- aggregatoR( data.bio )
#' df <- ambiguousSolver( data.agR )
#' df.bio <- asBiomonitor( df )
#' df.agR <- aggregatoR( df.bio )


ambiguousSolver <- function( x , method = "MCWP" ){

  classCheck( x )

  tree <- x[[ "Tree" ]]

  taxa <- as.character( tree$Taxa )
  taxa.pos <- which( names( tree ) == "Taxa" )

  # intialize the vector that will store the taxa already processed
  # this should speed up the computation time
  taxa.acc <- c()

  if( method == "RPKC" ){
    for( i in 1:length( taxa ) ){
      temp <- tree[  i , 1:(taxa.pos - 1) ]
      temp <- temp[ temp != "" ]
      temp <- temp[ -length( temp ) ]
      if( any( temp %in% taxa ) ){
        taxa.acc <- c( taxa.acc , temp[ temp %in% taxa ] )
      }
    }

    taxa.acc <- unique( taxa.acc )
    taxa.pc <- rep( TRUE, length( taxa ) )
    taxa.pc[ taxa %in% taxa.acc ] <- FALSE

    df <- x[["Taxa"]][ taxa.pc , ]
  }

  if( method == "MCWP"){
    for( i in 1:length( taxa ) ){taxa
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

    df <-  data.frame( Taxa = as.factor( taxa ) , tree[ , ( taxa.pos + 1 ):ncol( tree ) , drop = FALSE ] )

  }

  df
}
