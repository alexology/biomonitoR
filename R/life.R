#' life
#'
#' This function calculates LIFE index according to most recent version used in UK.
#' @param x results of function aggregatoR function
#' @param taxLev currently only the option "Family" is enabled
#' @param version possible choices are "extence" and "life_2017"
#' @param composite if TRUE composite families as listed in the details section are used
#' @param abucl abundance threshold. Default 0, 9, 99, 999, 9999
#' @keywords life
#' @details Lotic-invertebrate Index for Flow (LIFE) was originally proposed by Extence et al. (1999). biomonitoR implements the Extence et al. (1999) version called "extence" and the version currently used in UK called "life_2017".
#' If composite is set to T the following composite families are used for extence
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
#' while for "life_2017" the following are used:
#'
#'  \enumerate{
#'    \item Hydrophilidae (inc. Georissidae, Hlophoridae, Hydrochidae)
#'
#'  }
#' Scores used for life calculation can be explored with the function code{\link{showscores}}.
#'
#' @references Extence CA, Balbi DM, Chadd RP. 1999. River flow indexing using British benthic macroinvertebrates: a framework for setting hydroecological objectives. Regulated Rivers: Research and Management 15: 543-574.
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated LIFE scores and algorithms.
#' @importFrom stats aggregate reshape
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.life <- life(data.agR, taxLev = "Family", composite = FALSE)

life <- function(x, taxLev = "Family", version = "extence", composite = FALSE, abucl = c(0,9,99,999,9999)){

  # check if the object x is of class "biomonitoR"
  classCheck(x, group = "mi")

  if(taxLev != "Family"){
    stop("Species level LIFE not implemented yet")
  }


  if(taxLev == "Family" & version == "extence"){
    life_scores_use <- life_scores_fam
    fs_life_use <- fs_life
  }

  if(taxLev == "Family" & version == "life_2017"){
    life_scores_use <- life_scores_fam_2017
    fs_life_use <- fs_life
  }


  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  fam <- x[-numb, drop = FALSE]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites

  # composite default families

  # take into account composite families
  if(composite == TRUE & taxLev == "Family" & version == "extence"){
    fam <- checkBmwpFam(df=fam, famNames=life_fam_acc, stNames=st.names)
  }

  # take into account composite families
  if(composite == TRUE & taxLev == "Family" & version == "life_2017"){
    fam <- checkBmwpFam(df=fam, famNames=life_fam_acc_2017, stNames=st.names)
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
  temp <- fam.long[, 3, drop = FALSE]

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
  fam.long <- merge(fam.long, life_scores_use)
  names(fam.long)[5] <- "FS"
  fam.long <- merge(fam.long, fs_life_use)

  fam.sub <- fam.long[,c(4,9)]
  fam.life <- aggregate(. ~ Site, fam.sub, FUN = sum)
  fam.life$rich <- aggregate(. ~ Site, fam.sub, FUN = length)[,2]
  fam.life$score <- fam.life[,2] / fam.life[,3]
  res <- fam.life[, "score"]
  names(res) <- fam.life[, "Site"]
  res <- res[st.names]
  names(res) <- st.names
  return(res)
}
