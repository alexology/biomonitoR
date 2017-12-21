#' life
#'
#' This function calculates LIFE index according to Extence et al. (1999).
#' @param x results of function aggregatoR function
#' @param method possible choices
#' @keywords life
#' @details life currently calculates life index to family and species level, following the nomenclature used by Extence et al. (1999). Species level nomenclature of Extence et al. (1999) is outdated, an updated version will be released soon in biomonitoR.
#' @references Extence CA, Balbi DM, Chadd RP. 1999. River flow indexing using British benthic macroinvertebrates: a framework for setting hydroecological objectives. Regulated Rivers: Research and Management 15: 543–574.
#' @export
#' @seealso \code{\link{asBiomonitor}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' data.life <- life(data.agR, method = "Family")

life <- function(x, method = "Family"){

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  if(method == "Family"){
    life_scores <- life_scores_fam
  }

  if(method == "Species"){
    life_scores <- life_scores_spe
  }

  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  fam <- x[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites

  # take into account grouped families
  fam <- checkBmwpFam(df=fam, famNames=ukfam_acc, stNames=st.names)

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

  # transform row abundances to abunance classes
  names(temp) <- "ABU_NUM"
  temp[temp==0] <- 0
  temp[temp>=1&temp<=9] <- 1
  temp[temp>=10&temp<=99] <- 2
  temp[temp>=100&temp<=999] <- 3
  temp[temp>=1000&temp<=9999] <- 4
  temp[temp>=10000] <- 5

  fam.long <- data.frame(fam.long, temp)
  fam.long <- merge(fam.long, life_scores)
  names(fam.long)[5] <- "FS"
  fam.long <- merge(fam.long, fs_life)

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
