#' whpt
#'
#' This function calculates WHPT index according to most recent version used in UK.
#' @param x results of aggregatoR function
#' @param taxLev currently only the option "Family" is enabled
#' @param type possible choices are "aspt", "ntaxa", "bmwp"
#' @param composite if T composite families as listed in the details section are used
#' @param abucl abundance threshold. Default 0, 9, 99, 999, 9999.
#' @keywords whpt
#' @details WHPT is a revision of BMWP and it takes into account the abundances of organisms. The following aggregation is used if composit is set equal to T:
#'
#' \enumerate{
#'   \item Psychomyiidae (inc. Ecnomidae)
#'   \item Rhyacophilidae (inc. Glossomatidae)
#'   \item Ancylidae (inc. Acroloxidae)
#'   \item Gammaridae (inc. Crangonyctidae)
#'   \item Planariidae (inc. Dugesidae)
#'   \item Hydrobiidae (inc. Bithyniidae)
#' }
#'
#' Scores used for whpt calculation can be explored with the function code{\link{showscores}}.
#'
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated WHPT scores and algorithms.
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.whpt <- whpt(data.agR, taxLev = "Family", composite = F)

whpt <- function(x, taxLev = "Family", type = "aspt", composite = F, abucl = c(0,9,99,999,9999)){
  
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }
  
  if(type != "aspt" & type != "ntaxa" & type != "bmwp"){
    stop("Please provide a valide type: aspt, ntaxa or bmwp")
  }
  
  if(taxLev != "Family"){
    stop("Species level WHPT not implemented yet")
  }
  
  
  if(taxLev == "Family"){
    whpt_scores_use <- whpt_scores_fam
  }
  

  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  fam <- x[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
  
  # composite default families
  
  # take into account composite families
  if(composite == T & taxLev == "Family"){
    fam <- checkBmwpFam(df=fam, famNames=whpt_fam_acc, stNames=st.names)
  }
  
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
  E <- abucl[5]
  
  # transform row abundances to abunance classes
  names(temp) <- "ABU_NUM"
  temp[temp==A] <- 0
  temp[temp>=c(A+1) &temp<=B] <- 1
  temp[temp>=c(B+1)&temp<=C] <- 2
  temp[temp>=c(C+1)&temp<=D] <- 3
  temp[temp>=c(D+1)&temp<=E] <- 4
  temp[temp>=c(E+1)] <- 5
  
  fam.long <- data.frame(fam.long, temp)
  fam.long <- merge(fam.long, whpt_scores_use)

  fam.sub <- fam.long[,c(3,5)]
  fam.whpt <- aggregate(. ~ Site, fam.sub, FUN = sum)
  fam.whpt$rich <- aggregate(. ~ Site, fam.sub, FUN = length)[,2]
  fam.whpt$score <- fam.whpt[,2] / fam.whpt[,3]
  if(type == "aspt"){
    res <- fam.whpt[, "score"]
  }
  if(type == "ntaxa"){
    res <- fam.whpt[, "rich"]
  }
  if(type == "bmwp"){
    res <- fam.whpt[, "Score"]
  }  
  names(res) <- fam.whpt[, "Site"]
  res <- res[st.names]
  names(res) <- st.names
  return(res)
}
