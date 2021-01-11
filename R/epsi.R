#' e-psi
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates the Empyrically-weighted Proportion of Sediment-sensitive Invertebrates index (ePSI) according to the most recent version used in UK.
#' @param x results of the function `aggregate_taxa`.
#' @param method method `uk`.
#'  Users can provide their own data.frame (see examples) with a column called *Taxon* and the column of scores called *Scores*.
#' @param agg this option allows the composite family approach. It can be `FALSE`, `TRUE` or a `data.frame`.
#' If `FALSE` no aggregation will be performed, while if `TRUE` aggregation will be performed according to the rules described in Details.
#' A `data.frame` containing the aggregation rules can be provided by the user.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' `agg` cannot be `TRUE` when a `data.frame` is provided as method.
#' @param agg  optional, a data.frame provided by the user containing the specification on how to aggregate taxa. This data.frame needs a column called *Taxon*
#'  containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications. Used when users want to aggregate some taxonomic levels while using their own data.frame of scores (provided via `method`).
#' @param abucl Log abundance categories. Tresholds are set to 1, 9, 99, 999.
#' @param exceptions taxa that need to be exluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' `agg` cannot be `TRUE` when a data.frame is provided as `method`.
#' @param traceB if set to `TRUE` a list as specified below will be returned.
#' @keywords epsi
#' @details `epsi()` implementation take into account composite taxa as follow:
#' \enumerate{
#'  \item Tipulidae (inc. Limoniidae, Pediciidae & Cylindrotomidae)
#'  \item Siphlonuridae (inc. Ameletidae)
#'  \item Hydrophilidae (inc. Georissidae, Helophoridae & Hydrochidae)
#'  }
#'
#' The `epsi()` function automatically check for parent-child pairs in the scoring system, see the return section for a definition.
#' All the information used for `epsi()` calculation can be retrieved with the function \code{\link{show_scores}}.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Results of `epsi()`.
#'  \item `taxa_df` The data.frame used for the calculation containing the abundance of taxa receiving a score.
#'  \item `epsi_df` The data.frame used for the calculation containing scores and abundance classes for each site.
#'  \item `composite_taxa` Taxa aggregated following the aggregation rules of the `uk_agg` method or set in the `agg` option.
#'  \item `exceptions` A data.frame containing the containing changes made by excluding the taxa included in `exceptions`.
#'  \item `parent_child_pairs` For instance in Spanish `bmwp` both Ferrissia and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#'
#' @references Turley MD, Bilotta GS, Chadd RP, Extence CA, Brazier RE, Burnside NG, Pickwell AG. 2016. A sediment-specific family-level biomonitoring tool to identify the impacts of fine sediment in temperate rivers and streams. Ecological Indicators 70, 151-165.
#' @references Turley MD, Bilotta GS, Krueger T, Brazier RE, Extence CA. 2015. Developing an improved biomonitoring tool for fine sediment: combining expert knowledge and empirical data. Ecological indicators 54, 82-86.
#' @section Acknowledgements: We thank Carol Fitzpatrick, Richard Chadd, Judy England and Rachel Stubbington for providing us with the most updated ePSI scores and algorithms.
#' @importFrom stats aggregate reshape
#' @export
#' @seealso \code{\link{aggregate_taxa}}
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' epsi(data_agr)
#'
#' # provide your own score sistem. Scores and aggregation rules are for example purpose only.
#'
#' epsi_scores <- data.frame(
#'   Taxon = c("Ephemerellidae", "Leuctridae", "Chironomidae"),
#'   Scores = c(0.1, 0.5, 0.2)
#' )
#' epsi_acc <- data.frame(Taxon = "Ephemerellidae", Correct_Taxon = "Chironomidae")
#'
#' epsi(data_agr, method = epsi_scores, agg = epsi_acc, traceB = TRUE)
epsi <- function(x, method = "uk", agg = FALSE, abucl = c(1, 9, 99, 999), exceptions = NULL, traceB = FALSE) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # useful for transforming data to 0-1 later
  if (inherits(x, "bin")) {
    BIN <- TRUE
  } else {
    BIN <- FALSE
  }


  if (!identical(method, "uk") & !(is.data.frame(method))) (stop("Please provide a valid method"))
  if (!any(isFALSE(agg) | isTRUE(agg) | is.data.frame(agg))) stop("agg is not one of TRUE, FALSE or a custom data.frame")

  # Store tree for searching for inconsistencies
  Tree <- x[["Tree"]][, 1:10]

  numb <- c(which(names(x) == "Tree"), which(names(x) == "Taxa")) # position of the Tree and Taxa data.frame in the biomonitoR object that need to be removed

  # remove Tree and Taxa data.frame
  x <- x[-numb]
  st.names <- names(x[[1]][-1]) # names of the sampled sites
  # initialize the aggregation method
  z <- NULL


  # the following if statement is to allow the users to provide their own bmwp scores and aggregation rules.
  # y represents the method to be used
  if (is.data.frame(method) == TRUE) {
    if (isFALSE(agg)) {
      y <- method
    } else {
      y <- method
      z <- agg
    }
  } else {
    if (!(isTRUE(agg) | isFALSE(agg))) stop("When using the deafult method agg can only be TRUE or FALSE")

    # assign the default scores and aggregation rules as needed by the user

    if (identical(method, "uk")) {
      y <- epsi_scores_fam_uk

      if (isTRUE(agg)) {
        z <- epsi_acc_fam_uk
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

  if (BIN) {
    DF <- to_bin(DF)
  }

  if (traceB == "TRUE") {
    df2 <- DF
  }

  class.fun <- function(x) cut(x, breaks = c(abucl, 10^18), labels = 1:length(abucl), include.lowest = TRUE, right = TRUE)
  abu.class <- apply(apply(DF[, -1, drop = FALSE], 2, class.fun), 2, as.numeric)
  abu.class[is.na(abu.class)] <- 0
  DF <- data.frame(Taxon = DF[, 1], abu.class, check.names = FALSE)
  tot.mer <- merge(y, DF)

  EPSI <- data.frame(tot.mer[, "Scores"] * tot.mer[, st.names, drop = FALSE])
  epsi.sens <- apply(EPSI[tot.mer[, "Scores"] >= 0.5, , drop = FALSE], 2, sum)
  epsi.insens <- apply(EPSI, 2, sum)
  res <- epsi.sens / epsi.insens * 100
  names(res) <- st.names

  if (traceB == FALSE) {
    res
  } else {
    if (!exists("taxa.to.change")) {
      df3 <- NA
    } else {
      df3 <- taxa.to.change
    }
    if (exists("exce", inherits = FALSE)) {
      df4 <- exce
    } else {
      df4 <- "none"
    }
    if (exists("incon", inherits = FALSE)) {
      df5 <- incon
    } else {
      df5 <- "none"
    }
    res <- list(results = res, taxa_df = df2, epsi_df = tot.mer, composite_taxa = df3, exceptions = df4, parent_child_pairs = df5)
    res
  }
}
