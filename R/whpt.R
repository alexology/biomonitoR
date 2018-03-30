#' whpt
#'
#' This function calculates WHPT index according to most recent version used in UK.
#' @param x results of aggregatoR function
#' @param taxLev currently only the option "Family" is enabled
#' @param type presence only ("po") or abundance ("ab")
#' @param method possible choices are "aspt", "ntaxa", "bmwp"
#' @param composite if T composite families as listed in the details section are used
#' @param abucl abundance threshold. Default 0, 9, 99, 999.
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

whpt <- function(x, taxLev = "Family", type = "ab", method = "aspt", composite = F, abucl = c(0,9,99,999)){

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  if(type != "po" & type != "ab"){
    stop("Please provide a valide type: po or ab")
  }

  if(method != "aspt" & method != "ntaxa" & method != "bmwp"){
    stop("Please provide a valide metric: aspt, ntaxa or bmwp")
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



  # presence only

  if(type == "po"){
    whpt_scores_use <- whpt_scores_use[ which( whpt_scores_use$ABU_NUM == 1) , ] # keep the po scores only
    fam.long <- fam.long[ which( fam.long$Abu > 0) , ]
    fam.long <- merge(fam.long, whpt_scores_use)
  }

  if(type == "ab"){

    whpt_scores_use <- whpt_scores_use[ which( whpt_scores_use$ABU_NUM != 1) , ] # exclude the po scores

    # keep only numeric columns
    temp <- fam.long[, 3, drop = F]

    A <- abucl[1]
    B <- abucl[2]
    C <- abucl[3]
    D <- abucl[4]


    # transform row abundances to abunance classes
    names(temp) <- "ABU_NUM"
    temp[temp==A] <- 0
    temp[temp>=c(A+1) & temp<=B] <- 2
    temp[temp>=c(B+1) & temp<=C] <- 3
    temp[temp>=c(C+1) & temp<=D] <- 4
    temp[temp>=c(D+1)] <- 5


    fam.long <- data.frame(fam.long, temp)
    fam.long <- merge(fam.long, whpt_scores_use)
  }

  columns <- c("Site", "Score") # columns to retain for calculations
  fam.sub <- fam.long[ , columns]
  fam.whpt <- aggregate(. ~ Site, fam.sub, FUN = sum)
  fam.whpt$rich <- aggregate(. ~ Site, fam.sub, FUN = length)[,2]
  fam.whpt$score <- fam.whpt[,2] / fam.whpt[,3]
  if(method == "aspt"){
    res <- fam.whpt[, "score"]
  }
  if(method == "ntaxa"){
    res <- fam.whpt[, "rich"]
  }
  if(method == "bmwp"){
    res <- fam.whpt[, "Score"]
  }
  names(res) <- fam.whpt[, "Site"]
  res <- res[st.names]
  names(res) <- st.names
  return(res)
}
