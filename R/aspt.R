#' aspt
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates the *Average Score Per Taxon* index following Armitage et al. (1983), Davy-Bowker et al. (2007) and Alba-Tercedor & Sanchez-Ortega (1988) implementations.
#'
#' @param x Result of `aggregate_taxa()`.
#' @param method The implementation of BMWP needed to calculate ASPT. Possible choices are `a` (Armitage et al. 1983), `uk` (Davy-Bowker et al. 2010), `spa` (MAGRAMA 2011), `ita` (Buffagni et al . 2014).
#'  Users can provide their own data.frame (see examples) with a column called *Taxon* and the column of scores called *Scores*.
#' @param agg This option allows the composite family approach. It can be `FALSE`, `TRUE` or a `data.frame`.
#' If `FALSE` no aggregation will be performed, while if `TRUE` aggregation will be performed according to the rules described in Details.
#' A `data.frame` containing the aggregation rules can be provided by the user.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' `agg` cannot be `TRUE` when a `data.frame` is provided as method.
#' @param exceptions Taxa that need to be excluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' @param traceB If set to `TRUE` a list as specified below will be returned.
#'
#' @keywords aspt
#'
#' @references Armitage, P. D., Moss, D., Wright, J. F., & Furse, M. T. (1983). The performance of a new biological water quality score system based on macroinvertebrates over a wide range of unpolluted running-water sites. Water research, 17(3), 333-347.
#' @references Davy-Bowker J., Clarke R., Corbin T., Vincent H, Pretty J., Hawczak A., Blackburn J., Murphy J., Jones I., 2008. River Invertebrate Classification Tool. Final report. WFD72C. SNIFFER. 276 pp
#' @references MAGRAMA-Ministerio de Agricultura y medio Ambiente (2011) Protocolo de muestreo y laboratorio de fauna bentonica de invertebrados en rios vadeables. ML-Rv-I-2011, Cod, 23 pp.
#'
#' @details ASPT represents the average scores of the families that receive the score. Armitage scores are not reliable yet, since taxonomy has to be revised (e.g. Elminthidae are present instead of Elmidae). Davy-Bowker et al. (2010) and Buffagni et al. (2014) implementations take into account composite taxa as follow:
#' \enumerate{
#'   \item Psychomyiidae (inc. Ecnomidae)
#'   \item Rhyacophilidae (inc. Glossosomatidae)
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
#' Optional scores provided by the user data.frame needs to be formatted like following:
#' \tabular{lc}{
#' Taxon \tab Scores \cr
#' Aeshnidae \tab 8 \cr
#' Ancylidae \tab 6 \cr
#' Aphelocheiridae \tab 10 \cr
#' Asellidae \tab 3 \cr
#' Astacidae \tab 8 \cr
#' }
#'
#' Optional aggregation `data.frame` provided by the user needs to be formatted like following:
#' \tabular{ll}{
#'  Taxon \tab Correct_Taxon \cr
#'  Glossosomatidae \tab Rhyachopilidae \cr
#'  Apataniidae \tab Limnephilidae \cr
#'  Acroloxidae \tab Ancylidae \cr
#'  Crangonyctidae \tab Gammaridae \cr
#'  Niphargidae \tab Gammaridae \cr
#' }
#'
#'
#' The `aspt()` function automatically check for parent-child pairs in the scoring system, see the return section for a definition.
#' All the information used for `aspt()` calculation can be retrieved with the function \code{\link{show_scores}}.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of `aspt()`.
#'  \item `taxa_df` The data.frame used for the calculation containing the abundance of taxa receiving a score.
#'  \item `composite_taxa` Taxa aggregated following the aggregation rules when agg is not `NULL`.
#'  \item `exceptions` A data.frame containing the containing changes made by excluding the taxa included in `exceptions`.
#'  \item `parent_child_pairs` For instance in Spanish aspt both Ferrissia and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#'
#' @importFrom stats aggregate
#' @export
#' @seealso [aggregate_taxa] [bmwp]
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' aspt(data_agr)
#' aspt(data_agr, method = "spa")
aspt <- function(x, method = "ita", agg = FALSE, exceptions = NULL, traceB = FALSE) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # useful for transforming data to 0-1 later
  if (inherits(x, "bin")) {
    BIN <- TRUE
  } else {
    BIN <- FALSE
  }


  if (!any(identical(method, "ita"), identical(method, "spa"), identical(method, "a"), identical(method, "uk")) & !(is.data.frame(method))) (stop("Please provide a valid method"))
  if (!any(isFALSE(agg) | isTRUE(agg) | is.data.frame(agg))) stop("agg is not one of TRUE, FALSE or a custom data.frame")

  numb <- c(which(names(x) == "Tree"), which(names(x) == "Taxa")) # position of the Tree and Taxa data.frame in the biomonitoR object that need to be removed

  # Store tree for searching for inconsistencies
  Tree <- x[["Tree"]][, 1:10]

  # remove Tree and Taxa data.frame
  x <- x[-numb]
  st.names <- names(x[[1]][-1]) # names of the sampled sites
  # initialize the aggregation method
  z <- NULL

  # the following if statement is to allow the users to provide their own bmwp scores and aggregation rules.
  # y represents the method to be used
  if (is.data.frame(method)) {
    if (!(isFALSE(agg) | is.data.frame(agg))) {
      stop("When method is a data.frame agg needs to be FALSE or a data.frame containing the aggregation rules")
    }
    if (isFALSE(agg)) {
      y <- method
    } else {
      y <- method
      z <- agg
    }
  } else {
    if (!(isTRUE(agg) | isFALSE(agg))) stop("When using the deafult method agg can only be TRUE or FALSE")

    # assign the default scores and aggregation rules as needed by the user


    if (identical(method, "a")) (y <- aspt_scores_fam_armitage)

    if (identical(method, "ita")) {
      y <- aspt_scores_fam_ita

      if (isTRUE(agg)) {
        z <- aspt_acc_fam_ita
      }
    }

    if (identical(method, "spa")) {
      y <- aspt_scores_fam_spa
    }

    if (identical(method, "uk")) {
      y <- aspt_scores_fam_uk

      if (isTRUE(agg)) {
        z <- aspt_acc_fam_uk
      }
    }
  }

  # the calculation of the index in biomonitoR consists in rbind all the taxonomic levels
  # in the biomonitoR object that has been previously deprived of Taxa and Tree elements and then merge
  # it with the scores data.frame.
  # The first step is to change the column name of the first column of each data.frame to
  # an unique name

  for (i in 1:length(x)) {
    colnames(x[[i]])[1] <- "Taxon"
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
  if (!is.null(z)) {
    taxa.to.change <- as.character(DF$Taxon[DF$Taxon %in% z$Taxon])
    DF <- checkBmwpFam(DF = DF, famNames = z, stNames = st.names)
  } else {
    DF <- DF
  }

  DF <- manage_inconsistencies(DF = DF, Tree = Tree)
  if (!is.data.frame(DF)) {
    incon <- DF[[2]]
    DF <- DF[[1]]
  }

  # transform the data.frame from abundance to presence-absence if needed
  if (BIN) {
    DF <- to_bin(DF)
  }

  DF <- merge(y, DF)

  if (traceB == TRUE) {
    df1 <- DF
  }

  tot.mer <- data.frame(DF[, 1:2, drop = FALSE], (DF[, -c(1:2)] > 0) * 1)

  names(tot.mer)[-c(1, 2)] <- st.names # assign site names, the first wo columns are taxa and scores
  tot.st <- which(names(tot.mer) %in% st.names) # column numbers of the site columns
  ntaxa <- colSums(tot.mer[, -c(1:2), drop = F] == 1) # taxa richness, used as denominator in the aspt calculation
  tot.aspt <- apply(tot.mer$Scores * tot.mer[, tot.st, drop = F], 2, sum) / ntaxa # calculate the aspt as bmwp times the taxa richness


  # return the results
  if (!traceB) {
    tot.aspt
  } else {
    if (!exists("taxa.to.change", inherits = FALSE)) {
      df2 <- "none"
    } else {
      df2 <- taxa.to.change
    }
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


    list(results = tot.aspt, taxa_df = df1, composite_taxa = df2, exceptions = df3, parent_child_pairs = df4)
  }
}
