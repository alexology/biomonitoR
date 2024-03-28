#' Biologique Macrofitique en Riviere index
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("experimental") }
#'
#' This function calculates the IBMR index following Minciardi et al. (2009).
#'
#' @param x Result of `aggregate_taxa()`.
#' @param method The implementation of IBMR. the default method is `ita` (Minciardi et al., 2009).
#'  Users can provide their own data.frame (see examples) with a column called *Taxon* and the column of scores called *Scores*.
#' @param exceptions Taxa that need to be excluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' @param coverage_coef Are data already represented with coverage coefficients? Set to `TRUE` if this is not the case.
#' @param coverage_classes Thresholds to transform abundances into coverage classes.
#' @param traceB If set to `TRUE` a list as specified below will be returned.
#'
#' @keywords ibmr
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Result of `aspt()`.
#'  \item `taxa_df` The `data.frame` used for the calculation containing the abundance of taxa receiving a score.
#'  \item `exceptions` A `data.frame` containing the changes made by excluding the taxa included in `exceptions`.
#'  \item `parent_child_pairs` For instance in Spanish ASPT both *Ferrissia* and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#'
#'
#' @references Minciardi, M. R., Spada, C. D., Rossi, G. L., Angius, R., Orru, G., Mancini, L., Pace, G., Marcheggiani, S., Puccinelli, C. (2009). Metodo per la valutazione e la classificazione dei corsi d'acqua utilizzando la comunita delle macrofite acquatiche. ENEA Rapporto Tecnico RT/2009/23/ENEA.
#'
#' @details The function `ibmr()` calculates the index as follow:
#' \deqn{ibmr= \frac{\sum_{i=1}^{S} C{i} E{i} K{i}} {\sum_{i=1}^{S} E{i} K{i}}}
#'
#' where S is the number of taxa, i the ith taxon, Ci the sensitivity score, Ei the stenoecia coefficient
#' and Ki the coverage coefficient.
#'
#' Be careful when using data already transformed in coverage classes. The function `ibmr()` will sum
#' the classes if needed and this will result in wrong results. For example, if a taxon has only genus level
#' scores and two species of this genus are present, `ibmr()` will sum their coverage classes. This can result
#' in taxa having a coverage class higher than the highest coverage class allowed.
#'
#'
#' All the information used for `ibmr()` calculation can be retrieved with the function [show_scores].
#'
#' @importFrom stats aggregate
#' @export
#' @seealso [aggregate_taxa]
#' @examples
#' data(oglio)
#' oglio_asb <- as_biomonitor(oglio, group = "mf", FUN = bin)
#' oglio_agg <- aggregate_taxa(oglio_asb)
#' ibmr(oglio_agg)

ibmr <- function(x, method = "ita", exceptions = NULL, coverage_coef = FALSE, coverage_classes = c(0.1, 1, 10, 50), traceB = FALSE){

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # useful for transforming data to 0-1 later
  if (inherits(x, "bin")) {
    BIN <- TRUE
  } else {
    BIN <- FALSE
  }


  numb <- c(which(names(x) == "Tree"), which(names(x) == "Taxa")) # position of the Tree and Taxa data.frame in the biomonitoR object that need to be removed

  # Store tree for searching for inconsistencies
  Tree <- x[["Tree"]][, 1:10]

  # remove Tree and Taxa data.frame
  x <- x[-numb]
  st.names <- names(x[[1]][-1]) # names of the sampled sites


  for (i in 1:length(x)) {
    colnames(x[[i]])[1] <- "Taxon"
  }


  if (is.data.frame(method)){
    y <- method
  } else {
    y <- ibmr_scores
  }

  # rbind the data.frames representing a taxonomic level each
  # aggregate is not necessary here
  DF <- do.call("rbind", x)
  rownames(DF) <- NULL
  DF <- aggregate(. ~ Taxon, DF, sum)

  if (!is.null(exceptions)) {
    DF <- manage_exceptions(DF = DF, Tree = Tree, y = y, Taxon = exceptions)
    if (!is.data.frame(DF)) {
      exce <- DF[[2]]
      DF <- DF[[1]]
    }
  }


  # merge the new data.frame with the score data.frame and change
  # the names of the taxa according to the aggregation rules if needed
  DF <- merge(DF, y[, "Taxon", drop = FALSE])

  DF <- manage_inconsistencies(DF = DF, Tree = Tree)
  if (!is.data.frame(DF)) {
    incon <- DF[[2]]
    DF <- DF[[1]]
  }

  if(any(DF[, ! colnames(DF) %in% "Taxon", drop = FALSE] %in% 0:4)){
    message("Data different from 0, 1, 2, 3 and 4 detected. Probably they are not coverage coefficients.")
  }

  if(! coverage_coef){
    DF <- abundance_classes(DF, abundance_classes = coverage_classes)
  }

  DF_scores <- merge(y, DF[, "Taxon", drop = FALSE])

  multiplier <- DF_scores[, "Csi"] * DF_scores[, "Ei"]

  num <- sweep(DF[, ! colnames(DF) %in% "Taxon", drop = FALSE], 1, multiplier, "*")
  den <- sweep(DF[, ! colnames(DF) %in% "Taxon", drop = FALSE], 1, DF_scores[, "Ei"], "*")

  tot_ibmr <- apply(num, 2, sum) /apply(den, 2, sum)

  # return the results
  if (!traceB) {
    tot_ibmr
  } else {

    df1 <- data.frame(DF_scores, DF[, ! colnames(DF) %in% "Taxon", drop = FALSE])

    if (exists("exce", inherits = FALSE)) {
      df3 <- exce
    } else {
      df3 <- "none"
    }
    if (exists("incon", inherits = FALSE)) {
      df4 <- incon
    } else {
      df4 <- "none"
    }


    list(results = tot_ibmr, taxa_df = df1, exceptions = df3, parent_child_pairs = df4)
  }

}
