#' psi
#'
#' This function calculates Proportion of Sediment-sensitive Invertebrates index (PSI) according to the most recent version used in UK.
#' @param x results of aggregatoR function
#' @param taxLev currently only the option "Family" is allowed
#' @param abucl abundance threshold. Default 0, 9, 99, 999.
#' @keywords psi
#' @details 
#'
#' Scores used for whpt calculation can be explored with the function code{\link{showscores}}.
#' @references Extence CA, Chadd RP, England J, Dunbar MJ, Wood PJ, Taylor ED. 2013. The assessment of fine sediment accumulation in rivers using macro-invertebrate community response. River Research and Applications 29, 17–55.
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated WHPT scores and algorithms.
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.psi <- psi(data.agR, taxLev = "Family", composite = F)

psi <- function(x, taxLev = "Family", abucl = c(0,9,99,999)){
  
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  if(taxLev != "Family"){
    stop("Species level WHPT not implemented yet")
  }
  
  
  if(taxLev == "Family"){
    psi_scores_use <- psi_scores_fam
    fssr_psi_use <- fssr_psi_fam
  }
  
  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  fam <- x[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
  

  for(i in 1:length(fam)){
    colnames(fam[[i]])[1] <- "Taxon"
  }
  
  fam <- do.call( "rbind" , fam )
  rownames( fam ) <- NULL
  fam <- aggregate(. ~ Taxon, fam, FUN = sum)
  
  fam.long <- reshape(fam, direction="long", varying=list(names(fam)[-1]), v.names="Abu",
                      idvar="Taxon", times = names(fam)[-1], timevar = "Site")
  rownames(fam.long) <- NULL
 

  # keep only numeric columns
  temp <- fam.long[, 3, drop = F]
    
  A <- abucl[1]
  B <- abucl[2]
  C <- abucl[3]
  D <- abucl[4]
    
    
  # transform row abundances to abunance classes
  names(temp) <- "ABU_NUM"
  temp[temp==A] <- 0
  temp[temp>=c(A+1) & temp<=B] <- 1
  temp[temp>=c(B+1) & temp<=C] <- 2
  temp[temp>=c(C+1) & temp<=D] <- 3
  temp[temp>=c(D+1)] <- 4
    
    
  fam.long <- data.frame(fam.long, temp)
  fam.long <- merge(fam.long, psi_scores_use)
  names(fam.long)[5] <- "FSSR"
  fam.long <- merge(fam.long, fssr_psi_use)
  
  columns <- c("Site", "FSSR", "SCORE") # columns to retain for calculations
  fam.sub <- fam.long[ , columns]
  psi.ab <- aggregate(. ~ Site, fam.sub[which(fam.sub$FSSR == 1 | fam.sub$FSSR == 2) , ], FUN = sum)
  psi.abcd <- aggregate(. ~ Site, fam.sub, FUN = sum)
  
  # the next two lines ot overcome the problem of having 0 organisms belonign to 1 and 2 FSSR
  psi.mer <- merge(psi.ab, psi.abcd, by = "Site", all.x = T, all.y = T)
  psi.mer[is.na(psi.mer)] <- 0
  
  res <- psi.mer$SCORE.x/ psi.mer$SCORE.y
  names(res) <- psi.ab[, "Site"]
  res <- res[st.names]
  names(res) <- st.names
  return(res)
}
