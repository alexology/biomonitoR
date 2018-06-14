#' bqies
#'
#' This function calculates the Benthic Quality Index for Italian Lakes. Boggero et al. (2016), Italian classification method for macroinvertebrates in lakes. Method summary.
#' @param d results of function aggregatoR
#' @keywords aggregatoR
#' @details bqies is calculated....
#' @references Boggero A., Zaupa S., Cancellario T., Lencioni V., Marziali L., Rossaro B., (2016), Italian classification method for macroinvertebrates in lakes. Method summary. REPORT
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' bqies(data.agR)
#' bqies(data.agR, method = "i")

bqies <- function( d , method = "bqies") {
  
  # check if the object d is of class "biomonitoR"
  
  
  if (class(d) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  
  numb <- c(which(names(d)=="Tree"), which(names(d)=="Taxa")) # position of the Tree element in the list to remove
  x <- d[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  y <- bqies_scores
  st.names <- names(x[[1]][-1]) # names of sampled sites
  
  for(i in 1:length(x)){
    names(x[[i]])[1] <- "Taxon"
  }
  
  df <- do.call( "rbind" , x )
  rownames( df ) <- NULL
  df <- aggregate(. ~ Taxon, df, sum)
  d_merg <- merge( y , df )
  
  
  # check if merge results provided valid data.frame
  if( nrow(tot.mer) == 0 ){
    opt <- options( show.error.messages = T )
    on.exit( options( opt ) )
    return("No valid taxon provided")
  }
  else {
    # Take off the Score column in d_merge table
    d_merg.noWeight <- subset(d_merg, select = -c(Scores))
    # Calculate Log10+1 of all Taxa with Score
    d_merg.log <- log10(d_merg.noWeight[, -1]+1)
    # Sum all Log10+1 for column
    d_merg.sum <- colSums(d_merg.log)
    # Count the number of values != 0 in d_merg
    d_merg.no0 <- colSums(d_merg.log > 0)
    # Count Taxa not listed in Score table.
    spec_no_presence <- colSums(x[,-1] > 0) - d_merg.no0
    # percentage of Taxa density in Score table. If < of 75% the 
    # station point cannot be considered
    percent <- (1-(colSums(x[,-1])-colSums(d_merg.noWeight[,-1]))/colSums(x[,-1]))*100
    percent <- round(percent, digits = 1)
    # BQI
    bqi <- (colSums(d_merg.log*d_merg$Scores))/d_merg.sum
    bqi <- round(bqi, digits = 3)
    
    if(method == "bqi"){
      output <- list("BQI" = bqi, "N. Taxa without Score" = spec_no_presence, 
                     "% Taxa with Score"=percent)
    }
    if(method == "bqies"){
      # BQIES
      bqies <- bqi*log10(d_merg.no0 +  spec_no_presence +1)* colSums(x[,-1])/(colSums(x[,-1])+5)
      bqies <- round(bqies, digits = 3)
      output <- list("BQIES" = bqies,"N. Taxa without Score" = spec_no_presence, 
                     "% Taxa with Score"=percent)
    }
  }
  return( output )
}