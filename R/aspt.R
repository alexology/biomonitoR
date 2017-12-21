#' aspt
#'
#' This function calculates the Average Score Per Taxon following Armitage et al. (1983), Davy-Bowker et al. (2007) and Alba-Tercedor & Sanchez-Ortega (1988) implementations.
#' @param d results of function aggregatoR
#' @param method the formulation of BMWP needed to calculate ASPT. Possible choises are "a" (Armitage et al. 1983), "uk" (Davy-Bowker et al. 2010), "spa" (MAGRAMA 2011), "ita" (Buffagni et al . 2014). Methods "uk_agg"and "ita_agg" implement the composite family approach.
#' @keywords aggregatoR
#' @references Armitage, P. D., Moss, D., Wright, J. F., & Furse, M. T. (1983). The performance of a new biological water quality score system based on macroinvertebrates over a wide range of unpolluted running-water sites. Water research, 17(3), 333-347.
#' @references Davy-Bowker J., Clarke R., Corbin T., Vincent H, Pretty J., Hawczak A., Blackburn J., Murphy J., Jones I., 2008. River Invertebrate Classification Tool. Final report. WFD72C. SNIFFER. 276 pp
#' @references MAGRAMA-Ministerio de Agricultura y medio Ambiente (2011) Protocolo de muestreo y laboratorio de fauna bentónica de invertebrados en ríos vadeables. ML-Rv-I-2011, Cód, 23 pp.
#' @details ASPT represents the average scores of the families that receive the score. Armitage scores are not reliable yet, since taxonomy has to be revised (e.g. Elminthidae are present instead of Elmidae). Davy-Bowker implementation take into account composite taxa as follow:
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
#' @export
#' @seealso \code{\link{aggregatoR}}
#' @examples
#' data(macro_ex)
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#' aspt(data.agR)
#' aspt(data.agR, method = "i")

aspt <- function( d , method = "a") {

  # check if the object d is of class "biomonitoR"


  if (class(d) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }


  numb <- c(which(names(d)=="Tree"), which(names(d)=="Taxa")) # position of the Tree element in the list to remove
  x <- d[-numb, drop = F]
  # y is the reference data.set for bmwp calculation
  st.names <- names(x[[1]][-1]) # names of sampled sites

  if(method == "a") (y <- aspt_h)

  if(method == "ita") {y <- aspt_b
  z <- bfam_acc}

  if(method == "spa") { y <- aspt_i }

  if(method == "uk") {y <- aspt_uk
  z <- ukfam_acc}

  if(method == "ita_agg" || method == "uk_agg") (x <- checkBmwpFam(df=x, famNames=z, stNames=st.names))

  for(i in 1:length(x)){
    colnames(x[[i]])[1] <- "Taxon"
  }

  df <- do.call( "rbind" , x )
  rownames( df ) <- NULL
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
    names(tot.mer)[-c(1,2)] <- st.names
    tot.st <- which(names(tot.mer)%in%st.names)
    ntaxa <- colSums(tot.mer[, -c(1:2), drop = F] == 1)
    tot.aspt <- apply(tot.mer$Value*tot.mer[ , tot.st, drop=F], 2, sum)/ntaxa
  }
  return( round(tot.aspt, 3) )
}
