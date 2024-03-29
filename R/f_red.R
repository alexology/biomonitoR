#' Functional redundancy
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("maturing") }
#'
#' This function calculates the functional redundancy based on trait categories.
#'
#' @details
#' Functional redundancy (FR) is measured as the difference between taxonomic
#' diversity and functional diversity (de Bello et al., 2007). It relates positively
#' to ecosystem stability, resistance and resilience (Hooper et al. 2005; Guillemot
#' et al., 2011).
#'
#' \deqn{ FR = D - Q}
#'
#' The Gini-Simpson index is used to quantify Taxonomic Diversity (D, which ranges
#' from 0 to 1) where pi is the proportion of the abundance of taxa i in a
#' biological community.
#'
#' \deqn{D = 1 - \sum_{i=1}^{S} p_i^{2}}
#'
#' Rao quadratic entropy (Q; Rao, 1982) was used to estimate Functional Diversity
#' because it has been considered more appropriate than other indices (Botta-Dukat,
#' 2005; Ricotta, 2005). In this formula (Q) dij is the dissimilarity (ranging
#' from 0 to 1), between species i and j based on a set of specified functional
#' traits (i.e. effect traits, see below). This index is standardized by the maximum value to constrain the values
#' within the range of 0-1. Rao index is estimated using presence or abundance data
#' and the Euclidean transformed version of the traits-based Gower dissimilarity
#' matrix. For this, Gower's dissimilarity index (which ranges from 0 to 1) is used
#' because it can deal with traits of different nature and measuring scales
#' (continuous, nominal, binary, ordinal, etc.; see Podani 1999 for more information)
#'
#' \deqn{Q = \sum_{i=1}^{S} \sum_{j=1}^{S} d_{ij} \ p_i \ p_j}
#'
#' Given that the concept of FR was originally developed to represent the number
#' of taxa contributing similarly to an ecosystem function (Walker 1992; Lawton and
#' Brown 1993; Rosenfeld, 2002), FR and therefore functional diversity should be
#' calculated using only effect traits. Effect traits are those biological
#' features that directly influence a specific function of the ecosystem
#' (e.g. productivity, nutrient cycling). See Schmera et al. (2017) and Hevia et al (2017) for more
#' information about effect traits in aquatic invertebrate communities. Regarding
#' the interpretation of the results, when taxa within a community differ
#' completely in their functional traits, then Q = D and thus FR = 0. On the
#' other hand, when all taxa have identical functional traits, then Q = 0 and
#' FR = D, and when in addition the number of taxa is very large and they are
#' equally abundant, then D (and in this case FR) approaches 1 (Pillar et al., 2013).
#' Although the concept of FR could suggest that functionally similar species may
#' compensate for the loss or failure of others, there is evidence that ecosystems
#' need such redundancy to perform their functions efficiently and stably over time
#' (Rosenfeld 2002; Biggs et al. 2012). In fact, a decrease in FR could be dramatic
#' in non-redundant communities since the loss or replacement of one species could
#' lead to loss of unique traits or functions (Hooper et al. 2005), increasing
#' ecosystem vulnerability (Elmqvist et al. 2003).
#'
#' The `gower` distance refers to the mixed-variables coefficient of distance of Pavoine et al. (2009) as implemented in the `ade4` package.
#' This distance is meant to be used with fuzzy data.
#'
#' @param x Results of `aggregate_taxa()`.
#' @param trait_db A trait dataset. Can be a `data.frame` ot a `dist` object.
#' Taxonomic level of the traits dataset must match those of the taxonomic database.
#' No automatic check is done.
#' @param tax_lev Character string giving the taxonomic level used to retrieve
#' trait information. Possible levels are `Taxa`, `Species`, `Genus`,
#' `Family` as returned by `aggregate_taxa()`.
#' @param type The type of variables speciefied in `trait_db`.
#' Must be one of `F`, fuzzy, or `C`, continuous.
#' If more control is needed please consider to provide `trait_db` as a `dist` object.
#' It works only when `trait_db` is a `data.frame`, otherwise ingored.
#' @param traitSel Interactively select traits.
#' @param col_blocks A vector that contains the number of modalities for each trait.
#' Not needed when `euclidean` distance is used.
#' By default `biomonitoR` select the optimal number of dimensions with the quality of the functional space approach.
#' @param distance To be used to compute functional distances, `euclidean` or `gower`. Default to `gower`. See details.
#' @param zerodist_rm If `TRUE` aggregates taxa with the same traits.
#' @param correction Correction methods for negative eigenvalues, can be one of `none`, `lingoes` , `cailliez`, `sqrt` and `quasi`.
#' Ignored when type is set to `C`.
#' @param traceB if `TRUE` returns a list as specified in details.
#' @param set_param A list of parameters for fine tuning the calculations.
#' `tol` a tolerance threshold for zero, see the function `is.euclid`, `lingoes` and `cailliez` from the `ade4` for more details. Default to 1e-07.
#' If `cor.zero` is `TRUE`, zero distances are not modified. see the function `is.euclid`, `lingoes` and `cailliez` from the `ade4` for more details. Default to `TRUE`.
#'
#' @details Taxa without traits assigned in the trait database are removed from both the trait and abundance databases.
#'
#' @return a vector with fuzzy functional richness results.
#' \enumerate{
#'  \item `results` Results of `f_red()`.
#'  \item `traits` A `data.frame` containing the traits used for the calculations.
#'  \item `taxa` A `data.frame` conaining the taxa used for the calculations.
#'  \item `nbdim` Number of dimensions used after calculatin the quality of functional spaces according to Maire et al. (2015).
#'  \item `correction` The type of correction used.
#'  \item `NA_detection` A `data.frame` containing taxa on the first column and the corresponding trais with NAs on the second column.
#'  \item `duplicated_traits` If present, list the taxa with the same traits.
#' }
#'
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#' @importFrom stats complete.cases na.omit
#' @importFrom ade4 ktab.list.df dist.ktab prep.fuzzy divc quasieuclid is.euclid
#' @importFrom stats dist weighted.mean
#'
#' @examples
#' data(macro_ex)
#'
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#' # averaging
#' data_ts_av <- average_traits(data_ts)
#'
#' col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)
#'
#' f_red(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks)
#' f_red(data_agr,
#'   trait_db = data_ts_av, type = "F", col_blocks = col_blocks,
#'   correction = "cailliez"
#' )
#'
#' library(ade4)
#'
#' rownames(data_ts_av) <- data_ts_av$Taxa
#' traits_prep <- prep.fuzzy(data_ts_av[, -1], col.blocks = col_blocks)
#'
#' traits_dist <- ktab.list.df(list(traits_prep))
#' traits_dist <- dist.ktab(traits_dist, type = "F")
#'
#' f_red(data_agr, trait_db = traits_dist)
#' @seealso [aggregatoR]
#'
#' @references Biggs, R., Schluter, M., Biggs, D., Bohensky, E. L., BurnSilver, S.,
#'   Cundill, G., ... & Leitch, A. M. (2012). Toward principles for enhancing the
#'   resilience of ecosystem services. Annual Review of Environment and Resources,
#'   37, 421-448.
#' @references Botta-Dukat, Z. (2005). Rao's quadratic entropy as a measure of
#'   functional diversity based on multiple traits. Journal of Vegetation Science,
#'   16(5), 533-540.
#' @references de Bello, F., Leps, J., Lavorel, S., & Moretti, M. (2007).
#'   Importance of species abundance for assessment of trait composition:
#'   an example based on pollinator communities. Community Ecology, 8(2), 163-170.
#' @references Elmqvist, T., Folke, C., Nystrom, M., Peterson, G.,
#'   Bengtsson, J., Walker, B., & Norberg, J. (2003). Response diversity, ecosystem
#'   change, and resilience. Frontiers in Ecology and the Environment, 1(9), 488-494.
#' @references Guillemot, N., Kulbicki, M., Chabanet, P., & Vigliola, L. (2011).
#'   Functional redundancy patterns reveal non-random assembly rules in a
#'  species-rich marine assemblage. PLoS One, 6(10), e26735.
#' @references Hevia, V., Martin-Lopez, B., Palomo, S., Garcia-Llorente, M., de Bello, F.,
#'   & Gonzalez, J. A. (2017). Trait-based approaches to analyze links between the drivers
#'   of change and ecosystem services: Synthesizing existing evidence and future challenges.
#'   Ecology and evolution, 7(3), 831-844.
#' @references Hooper, D. U., Chapin, F. S., Ewel, J. J., Hector, A., Inchausti,
#'   P., Lavorel, S., et al. (2005). Effects of biodiversity on ecosystem
#'   functioning: a consensus of current knowledge. Ecological Monographs,
#'   75(1), 3-35.
#' @references Lawton, J.H. & Brown, V.K. (1993) Redundancy in ecosystems.
#'   Biodiversity and Ecosystem Function (eds E.-D. Schulze & H.A. Mooney),
#'   pp. 255-270. Springer-Verlag, Berlin.
#' @references Pavoine, S., Vallet, J., Dufour, A. B., Gachet, S., & Daniel, H. (2009).
#'  On the challenge of treating various types of variables: application for improving
#'   the measurement of functional diversity. Oikos, 118(3), 391-402.
#' @references Pillar, V. D., Blanco, C. C., Muller, S. C., Sosinski, E. E.,
#'   Joner, F., & Duarte, L. D. (2013). Functional redundancy and stability
#'   in plant communities. Journal of Vegetation Science, 24(5), 963-974.
#' @references Podani, J. (1999). Extending Gower's general coefficient of
#'   similarity to ordinal characters. Taxon, 331-340.
#' @references Rao, C. R. (1982). Diversity and dissimilarity coefficients:
#'   a unified approach. Theoretical population biology, 21(1), 24-43.
#' @references Ricotta, C. (2005). A note on functional diversity measures.
#'   Basic and Applied Ecology, 6(5), 479-486.
#' @references Rosenfeld, J. S. (2002). Functional redundancy in ecology
#'   and conservation. Oikos, 98(1), 156-162.
#' @references Schmera, D., Heino, J., Podani, J., Eros, T., & Doledec, S. (2017).
#'   Functional diversity: a review of methodology and current knowledge in
#'   freshwater macroinvertebrate research. Hydrobiologia, 787(1), 27-44.
#' @references Walker, B. H. (1992). Biodiversity and ecological redundancy.
#'   Conservation biology, 6(1), 18-23.
#'
#' @export

f_red <- function(x, trait_db = NULL, tax_lev = "Taxa", type = NULL, traitSel = FALSE, col_blocks = NULL, distance = "gower", zerodist_rm = FALSE, traceB = FALSE, correction = "none", set_param = NULL) {

  #  check if the object x is of class "biomonitoR"
  classCheck(x)

  # useful for transforming data to 0-1 later
  if (inherits(x, "bin")) {
    BIN <- TRUE
  } else {
    BIN <- FALSE
  }


  # list set_param for fine tuning of some functions
  if (is.null(set_param)) {
    set_param <- list(tol = 1e-07, cor.zero = TRUE)
  } else {
    set_param_def <- list(tol = 1e-07, cor.zero = TRUE)
    set_param_def[names(set_param)] <- set_param
    set_param <- set_param_def
  }

  if (is.null(trait_db)) stop("Please provide trait_db")

  if (!is.data.frame(trait_db) & !class(trait_db) %in% "dist") stop("trait_db must be a data.frame or a dist object")

  if (is.null(type) & is.data.frame(trait_db)) stop("Please specify a type when trait_db is a data.frame")

  if (!identical(type, "F") & !identical(type, "C") & is.data.frame(trait_db)) stop("type must be C or F when trait_db is a data.frame")

  if (identical(type, "C") & identical(distance, "gower")) (stop("Using gower distance when type is C is currently not allowed"))

  if (identical(type, "F") & identical(distance, "euclidean")) (warning("Are you sure to use euclidean distance when type is F?"))


  if (is.data.frame(trait_db)) {

    # trim and capitalise the first letter in the DB provided by the user
    trait_db[, "Taxa"] <- as.factor(sapply(trim(trait_db[, "Taxa"]), capWords, USE.NAMES = F))

    if (identical(type, "F")) {
      if (is.null(col_blocks)) (stop("Please provide col_blocks"))
      # check if the number of traits in trait_db equals the sum of col_blocks, otherwise stop
      if ((ncol(trait_db) - 1) != sum(col_blocks)) (stop("The number of traits in trait_db is not equal to the sum of col_blocks"))
    }
  }


  if (traitSel & is.data.frame(trait_db)) { # nocov start
    Index <- rep(1:length(col_blocks), col_blocks)
    rma <- select.list(names(trait_db[-which(names(trait_db) %in% "Taxa")]), title = "Traits selection", graphics = TRUE, multiple = T)
    # new col_blocks based on user trait selection, -1 because there is the column called Taxa
    col_blocks <- as.vector(table(Index[which(names(trait_db) %in% rma) - 1]))

    #  trait must have at least two modalities
    if (any(col_blocks < 2)) (stop("a trait must have at least two modalities"))

    trait_db <- trait_db %>%
      select(c("Taxa", rma))
    # trim and capitalise the column Taxa of the user' trait database
    trait_db$Taxa <- apply(as.data.frame(trim(trait_db$Taxa)), 1, capWords)
  } # nocov end


  st.names <- names(x[[1]][-1])

  DF <- x[[tax_lev]]
  names(DF)[1] <- "Taxon"

  taxa <- as.character(DF$Taxon)
  DF$Taxon <- as.character(DF$Taxon)


  if (is.data.frame(trait_db)) {
    trait_db$Taxa <- as.character(trait_db$Taxa)
    names(trait_db)[names(trait_db) %in% "Taxa"] <- "Taxon"

    # be sure that taxonomic and functional database have the same order and taxa
    DF <- merge(DF, trait_db[, "Taxon", drop = FALSE], by = "Taxon")


    # transform the data.frame from abundance to presence-absence if needed
    if (BIN) {
      DF <- to_bin(DF)
    }

    trait_db <- merge(trait_db, DF[, "Taxon", drop = FALSE], by = "Taxon")

    if (any(!DF$Taxon == trait_db$Taxon)) stop("Taxonomic and traits taxa does not match, ask the maintainer")

    # just to be sure we are doing the right things
    rownames(trait_db) <- trait_db$Taxon

    if (identical(type, "F")) (tr_prep <- prep.fuzzy(trait_db[, -1], col.blocks = col_blocks))
    if (identical(type, "B")) (tr_prep <- prep.binary(trait_db[, -1], col.blocks = col_blocks))
    if (identical(type, "C")) (tr_prep <- trait_db[, -1])

    rownames(tr_prep) <- trait_db$Taxon


    # computing functional dissimilarity between species given their traits values
    if (identical(distance, "gower")) {
      mat_dissim <- ktab.list.df(list(tr_prep))
      mat_dissim <- dist.ktab(mat_dissim, type = "F")
    }

    if (identical(distance, "euclidean")) mat_dissim <- dist(scale(tr_prep)) # scaling if continuous traits
  }

  if (class(trait_db) %in% "dist") {

    # keep only common taxa between taxonomic and traits database
    # to do this the dist object is transformed into matrix
    trait_db <- as.matrix(trait_db)
    DF <- merge(DF, data.frame(Taxon = rownames(trait_db)), by = "Taxon")

    # transform the data.frame from abundance to presence-absence if needed
    if (BIN) {
      DF <- to_bin(DF)
    }

    trait_db <- trait_db[rownames(trait_db) %in% DF$Taxon, ]
    trait_db <- trait_db[, colnames(trait_db) %in% DF$Taxon, drop = FALSE]
    trait_db <- trait_db[match(DF$Taxon, rownames(trait_db)), ]
    trait_db <- trait_db[, match(DF$Taxon, colnames(trait_db)), drop = FALSE]

    # check if names are in the same order, both on rows and columns
    if (any(!DF$Taxon == rownames(trait_db))) stop("Taxonomic and traits taxa does not match, ask the maintainer")
    if (any(!DF$Taxon == colnames(trait_db))) stop("Taxonomic and traits taxa does not match, ask the maintainer")

    mat_dissim <- as.dist(trait_db)
  }

  if (zerodist_rm & any(mat_dissim < set_param$tol)) {
    zero_corr <- zero_dist_traits(x = DF, mat_dissim = mat_dissim, BIN = BIN)
    DF <- zero_corr[[1]]
    mat_dissim <- zero_corr[[2]]
    df1 <- zero_corr[[3]]
  }


  if (identical(distance, "gower")) {
    if (identical(correction, "cailliez")) mat_dissim <- suppressWarnings(cailliez(mat_dissim, tol = set_param$tol, cor.zero = set_param$cor.zero))
    if (identical(correction, "lingoes")) mat_dissim <- suppressWarnings(lingoes(mat_dissim, tol = set_param$tol, cor.zero = set_param$cor.zero))
    if (identical(correction, "sqrt")) mat_dissim <- suppressWarnings(sqrt(mat_dissim))
    if (identical(correction, "quasi")) mat_dissim <- suppressWarnings(quasieuclid(mat_dissim))
  }

  suppressWarnings(euclid.dist.mat <- is.euclid(mat_dissim, tol = set_param$tol))

  if (!euclid.dist.mat) {
    stop("Non euclidean trait distance. Euclidean property is needed. Please use the correction options
             otherwise consider to remove taxa with the same traits.")
  }


  if (any(mat_dissim < set_param$tol)) {
    MES <- "At least a pair of species has the same traits. Depending on your needs, this could be an issue."
    message(MES)
  } else {
    MES <- "no taxa with the same traits"
  }


  rownames(DF) <- DF[, "Taxon"]

  tax_sim <- suppressWarnings(divc(DF[, -1])$diversity)
  raoQ <- suppressWarnings(divc(DF[, -1], mat_dissim, scale = T)$diversity)


  FRed <- tax_sim - raoQ
  FRed[FRed < 0] <- 0

  res <- data.frame(GS_rich = tax_sim, raoQ = raoQ, fred = FRed)
  rownames(res) <- st.names

  if (!traceB) {
    return(res)
  }

  if (traceB) {
    # chech for NA, it could happen that a trait is filled with NAs
    # but this can be done only when trait_db is a data.frame
    if (is.data.frame(trait_db)) {
      if (any(is.na(tr_prep))) {
        tax.na <- as.data.frame(which(is.na(tr_prep), arr.ind = TRUE))
        tax.na[, 1] <- trait_db[tax.na[, 1], 1]
        tax.na[, 2] <- colnames(trait_db[, -1])[tax.na[, 2]]
        colnames(tax.na) <- c("Taxa", "Traits")
        tax.na <- tax.na[order(tax.na[, 1]), ]
        rownames(tax.na) <- NULL
      } else {
        tax.na <- "No NAs detected"
      }
    } else {
      tax.na <- "NAs cannot be detected when trait_db is a dist object"
    }

    # prepare traits to be returned
    if (!is.data.frame(trait_db)) {
      # returns the distance matrix used for the calculation as a dist object
      trait_db <- mat_dissim
    }

    # prepare traits to be returned
    if (is.data.frame(trait_db)) {
      # returns the distance matrix used for the calculation as a dist object
      rownames(trait_db) <- NULL
    }


    if (exists("df1", inherits = FALSE)) {
      df1 <- df1
    } else {
      df1 <- MES
    }


    rownames(DF) <- NULL

    res.list <- list(res, trait_db, DF, correction = correction, tax.na, df1)
    names(res.list) <- c("results", "traits", "taxa", "correction", "NA_detection", "duplicated_traits")
    return(res.list)
  }
}
