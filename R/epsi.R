#' e-psi
#'
#' This function calculates Empyrically-weighted Proportion of Sediment-sensitive Invertebrates index (ePSI) according to the most recent version used in UK.
#' @param x results of aggregatoR function
#' @param taxLev currently only the option "Family" is allowed
#' @param abucl Log abundance categories. Treshold are set to 0, 9, 99, 999.
#' @param composite if T composite families as listed in the details section are used.
#' @keywords epsi
#' @details ePSI implementation take into account composite taxa as follow:
#' \enumerate{
#'  \item Tipulidae (inc. Limoniidae, Pediciidae & Cylindrotomidae)
#'  \item Siphlonuridae (inc. Ameletidae)
#'  \item Hydrophilidae (inc. Georissidae, Helophoridae & Hydrochidae)
#'  }
#' Ancylus (0.51) has a different score compared to the other Planorbidae (0). Take it into account when interpreting the results.
#' Scores used for epsi calculation can be explored with the function code{\link{showscores}}.
#' @references Turley MD, Bilotta GS, Chadd RP, Extence CA, Brazier RE, Burnside NG, Pickwell AG. 2016. A sediment-specific family-level biomonitoring tool to identify the impacts of fine sediment in temperate rivers and streams. Ecological Indicators 70, 151-165.
#' @references Turley MD, Bilotta GS, Krueger T, Brazier RE, Extence CA. 2015. Developing an improved biomonitoring tool for fine sediment: combining expert knowledge and empirical data. Ecological indicators 54, 82-86.
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated ePSI scores and algorithms.
#' @importFrom stats aggregate reshape
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.psi <- epsi(data.agR, taxLev = "Family", composite = F)

epsi <- function(x, taxLev = "Family", abucl = c(0,9,99,999), composite = F){

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  if(taxLev != "Family"){
    stop("Currently only family level is implemented!")
  }


  if(taxLev == "Family"){
    epsi_scores_use <- epsi_scores_fam
    epsi_fam_acc <- epsi_fam_acc
  }

  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  fam <- x[-numb, drop = F]

  st.names <- names(x[[1]][-1]) # names of sampled sites

    # take into account composite families
  if(composite == T & taxLev == "Family"){
    fam <- checkBmwpFam(df=fam, famNames=epsi_fam_acc, stNames=st.names)
  }

  for(i in 1:length(fam)){
    colnames(fam[[i]])[1] <- "Taxon"
  }

  fam <- do.call( "rbind" , fam )
  rownames( fam ) <- NULL
  fam <- aggregate(. ~ Taxon, fam, FUN = sum)

  # take into account the Ancylus/Planorbidae problem
  ancy <- fam[ which(fam$Taxon == "Ancylus") , ] # subset Ancylus

  if(nrow(fam[ which(fam$Taxon == "Ancylidae") , ]) > 0){ fam <- fam} else { # this row because if the user provide a custom database with Ancylidae the abundance of Ancylus must not be subtract from the abundance of Planorbidae

  if(nrow(ancy) > 0){
    taxon.col <- which(names(fam) %in% c("Taxon")) # Identifying the Taxon column, should be the first
    plan.row <- rownames(fam[ which(fam$Taxon == "Planorbidae") , ]) # identify the Planorbidae row
    fam[plan.row,-1] <- fam[ plan.row , -taxon.col] - ancy[ , -taxon.col] # -1 for not considering the column Taxon
    }
    else (fam <- fam)
  }

  fam.long <- reshape(fam, direction="long", varying=list(names(fam)[-1]), v.names="Abu",
                      idvar="Taxon", times = names(fam)[-1], timevar = "Site")
  rownames(fam.long) <- NULL


  # keep only numeric columns
  temp <- fam.long[, 3, drop = F]

  A <- abucl[1]
  B <- abucl[2]
  C <- abucl[3]
  D <- abucl[4]


  # transform row abundances to abunance categories
  names(temp) <- "ABU_NUM"
  temp[temp==A] <- 0
  temp[temp>=c(A+1) & temp<=B] <- 1
  temp[temp>=c(B+1) & temp<=C] <- 2
  temp[temp>=c(C+1) & temp<=D] <- 3
  temp[temp>=c(D+1)] <- 4


  fam.long <- data.frame(fam.long, temp)
  fam.long <- merge(fam.long, epsi_scores_use)
  fam.long$EPSI <- fam.long$ABU_NUM * fam.long$Scores
  epsi.sensitive <- aggregate(EPSI ~ Site, data = fam.long[which(fam.long$Scores >= 0.5),], FUN = sum)
  epsi.insensitive <- aggregate(EPSI ~ Site, data = fam.long, FUN = sum)

  # the next two lines ot overcome the problem of having 0 organisms belonging to 1 and 2 FSSR
  epsi.mer <- merge(epsi.sensitive, epsi.insensitive, by = "Site", all.x = T, all.y = T)
  epsi.mer[is.na(epsi.mer)] <- 0

  # x stands for sensitive epsi and y for insensitive epsi. This results from the merge function above
  res <- epsi.mer$EPSI.x/epsi.mer$EPSI.y*100
  names(res) <- epsi.mer[, "Site"]
  res <- res[st.names]
  names(res) <- st.names
  return(res)
}
