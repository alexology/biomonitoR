#' indNumb
#'
#' Function for calculating the number of Taxa according to taxalist of desired index.
#' @param d results of function aggregatoR
#' @param method the only option is "b" (see details).
#' @keywords indNumb
#' @details Only implemented for the Italian version of the STAR_ICMi
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{aggregatoR}}, \code{\link{speNumb}}, \code{\link{genNumb}}, \code{\link{famNumb}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' indNumb(data.agR)


indNumb <- function( d , method = "a") {

  # check if the object d is of class "biomonitoR"


  if (class(d) != "biomonitoR") {
    stop("The object is not of class biomonitoR")
  }


  numb <- c(which(names(d)=="Tree"), which(names(d)=="Taxa")) # position of the Tree element in the list to remove
  x <- d[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
  if(method == "b") {y <- aspt_b_fam}

  for(i in 1:length(x)){
    colnames(x[[i]])[1] <- "Taxon"
  }

  df <- do.call( "rbind" , x )
  rownames( df ) <- NULL
  df <- aggregate(. ~ Taxon, df, sum)
  tot.mer <- merge( y , df )

  # Column 1 represents Taxa names and it need to be excluded from the calculations
  ntax <- apply(tot.mer[, -1, drop = F], 2, FUN = function(x) { length( x[x>0] ) } )
  names(ntax) <- st.names
  return(ntax)

}
