#' bmwp
#'
#' This function calculates the Biological Monitoring Working Party following Armitage et al. (1983), Davy-Bowker et al. (2007) and Alba-Tercedor & Sanchez-Ortega (1988) implementations.
#'
#' @param x results of aggregatoR function.
#' @param method the implementation of BMWP needed to calculate ASPT. Possible choices are "a" (Armitage et al. 1983), "uk" (Davy-Bowker et al. 2010), "spa" (MAGRAMA 2011), "ita" (Buffagni et al . 2014). Methods "uk_agg"and "ita_agg" implement the composite family approach.
#'  Users can provide their own a data.frame (see examples) with a column called *Taxon* and the column of scores called *Value*.
#' @param agg a data.frame containing the specification on how to aggregate taxa. This data.frame needs a column called *Taxon*
#'  containing the taxon to aggregate and a column called *Taxon_changed* with the aggregation specifications.
#' @keywords aggregatoR
#' @references Armitage, P. D., Moss, D., Wright, J. F., & Furse, M. T. (1983). The performance of a new biological water quality score system based on macroinvertebrates over a wide range of unpolluted running-water sites. Water research, 17(3), 333-347.
#' @references Davy-Bowker J., Clarke R., Corbin T., Vincent H, Pretty J., Hawczak A., Blackburn J., Murphy J., Jones I., 2008. River Invertebrate Classification Tool. Final report. WFD72C. SNIFFER. 276 pp
#' @references MAGRAMA-Ministerio de Agricultura y medio Ambiente (2011) Protocolo de muestreo y laboratorio de fauna bentonica de invertebrados en rios vadeables. ML-Rv-I-2011, Cod, 23 pp.
#' @details BMWP is calculated as the sum of scores of the sensitive taxa present in a given sample. Armitage scores are not reliable yet, since taxonomy has to be revised (e.g. Elminthidae are present instead of Elmidae). Davy-Bowker et al. (2010) and Buffagni et al. (2014) implementations take into account composite taxa as follow:
#' \enumerate{
#'   \item Psychomyiidae (inc. Ecnomidae)
#'   \item Rhyachopilidae (inc. Glossomatidae)
#'   \item Limnephilidae (inc. Apatanidae)
#'   \item Ancylidae (inc. Acroloxidae)
#'   \item Gammaridae (inc. Crangonyctidae & Niphargidae)
#'   \item Hydrophilidae (inc. Hydraenidae, Helophoridae)
#'   \item Tipulidae (inc. Limoniidae, Pediciidae & Cylindrotomidae)
#'   \item Planariidae (inc. Dugesidae)
#'   \item Hydrobiidae (inc. Bithyniidae)
#'   \item Oligochaeta (all the families)
#' }
#'
#' User provided scores data.frame needs to be formatted like following:
#' \tabular{lc}{
#' Taxon \tab value \cr
#' Aeshnidae \tab 8 \cr
#' Ancylidae \tab 6 \cr
#' Aphelocheiridae \tab 10 \cr
#' Asellidae \tab 3 \cr
#' Astacidae \tab 8 \cr
#' }
#'
#' User provide aggregation data.frame needs to be formatted like following:
#' \tabular{ll}{
#'  Taxon \tab Taxon_changed \cr
#'  Glossomatidae \tab Rhyachopilidae \cr
#'  Apatanidae \tab Limnephilidae \cr
#'  Acroloxidae \tab Ancylidae \cr
#'  Crangonyctidae \tab Gammaridae \cr
#'  Niphargidae \tab Gammaridae \cr
#' }
#'
#' @note Carefully check if your taxa list contains Ancylidae. Ancylidae is considered
#' as family in italian and uk ASPT, while it is not considered a valid family by [freshwaterecology.info](https://www.freshwaterecology.info/).
#'
#' @importFrom stats aggregate
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' bmwp(data.agR)
#' bmwp(data.agR, method = "spa")


bmwp <- function( x , method = "ita", agg = NULL) {

  # check if the object x is of class "biomonitoR"
  classCheck(x, group = "mi")

  numb <- c(which(names(x)=="Tree"), which(names(x)=="Taxa")) # position of the Tree element in the list to remove
  x <- x[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites
  z <- NULL

  # the following if statement is to allow the users to provide their own bmwp scores and aggregation.
  if( is.data.frame(method ) == TRUE){
    if( is.null( agg ) == TRUE ){
      y <- method
    } else {
      y <- method
      z <- agg
    }
  }

  else{
    if(method == "a") (y <- aspt_h)

    if(method == "ita") (y <- aspt_b)

    if( method == "ita_agg" ) {y <- aspt_b
    z <- bfam_acc}

    if(method == "spa") { y <- aspt_i }

    if( method == "uk" ) ( y <- aspt_uk )

    if( method == "uk_agg") {y <- aspt_uk
    z <- bfam_acc_uk}
  }

  for(i in 1:length(x)){
    colnames(x[[i]])[1] <- "Taxon"
  }

  df <- do.call( "rbind" , x )
  rownames( df ) <- NULL
  df <- aggregate(. ~ Taxon, df, sum)
  df <- merge( df, y[ , "Taxon", drop = FALSE])
  if(is.null(z) == FALSE){
    df <- checkBmwpFam( df = df , famNames = z , stNames = st.names )
  } else {
    df <- df
  }
  df <- data.frame( df[ , 1 , drop =F ], (df[ , -1 ] > 0 ) * 1 )
  tot.mer <- merge( y , df )

  # check if merge results provided valid data.frame
  if( nrow(tot.mer) == 0 ){
    opt <- options( show.error.messages = T )
    on.exit( options( opt ) )
    return("No valid taxon provided")
  }
  else {
    names(tot.mer)[-c(1,2)] <- st.names
    tot.st <- which(names(tot.mer)%in%st.names)
    tot.bmwp <- apply(tot.mer$Value*tot.mer[ , tot.st, drop=F], 2, sum)
  }
  return( tot.bmwp )
}
