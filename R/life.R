#' life
#'
#' This function calculates LIFE index according to Davy-Bowker et al. (2010).
#' @param x results of function aggregatoR function
#' @param taxLev possible choices are "Family" and "Species"
#' @param composite if T composite families as listed in the details section are used
#' @param abucl abundance threshold. Default 0, 9, 99, 999, 9999
#' @keywords life
#' @details Lotic-invertebrate Index for Flow (LIFE) was originally proposed by Extence et al. (1999). biomonitoR implements the version proposed by Davy-Bowker et al. (2010) that has an updated taxonomic list compared to those in Extence et a. (1999). The life function aggregate the following families by default:
#' \enumerate{
#'   \item Limnephilidae (inc. Apatanidae)
#'   \item Gammaridae (inc. Niphargidae)
#'   \item Hydrophilidae (inc. Helophoridae, Georissidae & Hydrochidae)
#'   \item Tipulidae (inc. Limoniidae, Pediciidae & Cylindrotomidae)
#'   \item Siphlonuridae (inc. Ameletidae)
#' }
#'
#' If composite is set to T the following composite families are used
#'
#' \enumerate{
#'   \item Psychomyiidae (inc. Ecnomidae)
#'   \item Rhyachopilidae (inc. Glossomatidae)
#'   \item Ancylidae (inc. Acroloxidae)
#'   \item Gammaridae (inc. Crangonyctidae)
#'   \item Hydrophilidae (inc. Hydraenidae)
#'   \item Planariidae (inc. Dugesidae)
#'   \item Hydrobiidae (inc. Bithyniidae)
#'   \item Dytiscidae (inc. Noteridae)
#'   \item Planariidae (inc. Dugesiidae)
#' }
#'
#' Scores used for life calculation can be explored with the function code{\link{showscores}}.
#'
#' @references Extence CA, Balbi DM, Chadd RP. 1999. River flow indexing using British benthic macroinvertebrates: a framework for setting hydroecological objectives. Regulated Rivers: Research and Management 15: 543–574.
#' @references Davy-Bowker J, Arnott S, Close R, Dobson M, Dunbar M, Jofre G, Morton D, Murphy J, Wareham W, Smith S, Gordon V. 2010. SNIFFER WFD 100: Further development of River Invertebrate Classification Tool. Final Report.
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.life <- life(data.agR, taxLev = "Family", composite = F)

life <- function(x, taxLev = "Family", composite = F, abucl = c(0,9,99,999,9999)){

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  if(taxLev == "Family"){
    life_scores_use <- life_scores_fam
    fs_life_use <- fs_life
  }

  if(taxLev == "Species"){
    life_scores_use <- life_scores_spe
    fs_life_use <- fs_life
  }


  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  fam <- x[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites

  # composite default families

  fam <- checkBmwpFam(df=fam, famNames=life_acc_default, stNames=st.names)

  # take into account composite families
  if(composite == T){
    fam <- checkBmwpFam(df=fam, famNames=life_fam_acc, stNames=st.names)
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
  fam.long <- merge(fam.long, life_scores_use)
  names(fam.long)[5] <- "FS"
  fam.long <- merge(fam.long, fs_life_use)

  fam.sub <- fam.long[,c(4,9)]
  fam.life <- aggregate(. ~ Site, fam.sub, FUN = sum)
  fam.life$rich <- aggregate(. ~ Site, fam.sub, FUN = length)[,2]
  fam.life$score <- fam.life[,2] / fam.life[,3]
  res <- round(fam.life[, "score"], 3)
  names(res) <- fam.life[, "Site"]
  res <- res[st.names]
  names(res) <- st.names
  return(res)
}
