#' bmwp
#'
#' Functions for calculating BMWP and ASPT
#' @param x results of function aggregatoR
#' @param method a,b or i. See details.
#' @keywords aggregatoR
#' @details
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' bmwp(data.agR)



bmwp <- function( d , method = "a") {

  # check if the object d is of class "biomonitoR"


    if (class(d) != "biomonitoR") {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      return("Object x is not an object of class biomonitoR")
    }


  x <- d[-numb]

  
  numb <- c(which(names(d)=="Tree"), which(names(d)=="Taxa")) # position of the Tree element in the list to remove
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
  if(method == "a") (y <- aspt_h)
  if(method == "b") {y <- aspt_b
  z <- bfam_acc}
  if(method == "i") {y <- aspt_i
  z <- ifam_acc
  # solving the Planorbidae/Ancylidae problem
  # Ancylidae genus http://eol.org/pages/2402/names
  ancylidae.gen <- c("Ancylus", "Ferrissia", "Gundlachia", "Hebetoncylus", "Laevapex", "Rhodacmaea", "Rhodacme")
  temp.tree <- d[["Tree"]]
  ancylidae.sub <- temp.tree[which(temp.tree$Genus %in% ancylidae.gen), st.names, drop=F]
  ancylidae.abu <- apply(ancylidae.sub, 2, sum)
  ancylidae.pa <- ancylidae.abu
  ancylidae.pa[ which( ancylidae.pa >0 ) ] <-1
  }
  
  if(method == "b" || method == "i") (x <- checkBmwpFam(df=x, famNames=z, stNames=st.names))
  
  for(i in 1:length(x)){
    colnames(x[[i]])[1] <- "Taxon"
  }

  df <- do.call( "rbind" , x )
  rownames( df ) <- NULL
  
  # # solving the Planorbidae/Ancylidae problem
  if(method == "i"){
    planorbidae.row <- which(df$Taxon == "Planorbidae")
    ancylidae.row <- which(df$Taxon == "Ancylidae")
    if( length(planorbidae.row) != 0){
      df[planorbidae.row ,-1] <- df[planorbidae.row ,-1] - ancylidae.pa
    }
    if( length(ancylidae.row) != 0){
      temp.anc <- as.numeric(df[ancylidae.row ,-1]) + ancylidae.pa
      temp.anc[temp.anc > 0] <- 1
      df[ancylidae.row ,-1] <- temp.anc
    } else {
      levels(df$Taxon) <- c(levels(df$Taxon), "Ancylidae")
      temp.anc <- data.frame(Taxon = "Ancylidae", as.data.frame(t(ancylidae.pa)))
      df <- rbind(df, temp.anc)
    }
  }
  
  df <- aggregate(. ~ Taxon, df, sum)
  df <- data.frame( df[ , 1 , drop =F ], (df[ , -1 ] > 0 ) * 1 )
  tot.mer <- merge( y , df )

  # check if merge results provided valid data.frame
  if( nrow(tot.mer) == 0 ){
    opt <- options( show.error.messages = T )
    on.exit( options( opt ) )
    return("No valid taxon provided")
  }
  else {
    tot.st <- which(names(tot.mer)%in%st.names)
    tot.bmwp <- apply(tot.mer$Value*tot.mer[ , tot.st, drop=F], 2, sum)
  }
  return( tot.bmwp )
}
