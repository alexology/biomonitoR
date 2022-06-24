#' @name fuzzy_trait_ratio
#' @title Indices from fuzzy coded data
#'
#' @description This function calculates indices starting from fuzzy coded data.
#' @param x Result of `aggregate_taxa()`.
#' @param tax_lev Taxonomic level at which the calculation has to be performed.
#' @param trait_db A trait `data.frame` with a column Taxa and the other columns containing the traits. Please check `assign_traits()` for building the trait dataset.
#' @param type Presence only `po` or abundance `ab`.
#' @param numerator Names of the columns to set as the numerator.
#' @param denominator Names of the columns to set as the denominator.
#' @param trans The function used to transform the abundances, by default `log1p`.
#' @details This function best performs with standardized fuzzy coded data that sum up to 1 for each taxon.
#' For each taxon, it works by i) summing the scores at the numerator, ii) summing the scores at the denominator,
#' iii) multiplying these sums for the abundance of the target taxon. This step is performed for all the taxa in a community.
#' The last step consists in summing up the scores of all the taxa for both the numerator and denominator and
#' in taking the ratio between the numerator and the denominator.
#' An example of this approach is the Flow_T index developed by Laini et al. (2022).
#' @seealso [as_biomonitor] [assign_traits]
#' @export
#' @references Laini, A., Burgazzi, G., Chadd, R., England, J., Tziortzis, I., Ventrucci, M., Vezza, P., Viaroli, P., Wood, P.J. & Guareschi, S. (2022). Using invertebrate functional traits to improve flow variability assessment within European rivers. Science of The Total Environment, 832, 155047.
#' @examples
#'
#' data(mi_prin)
#' data_bio <- as_biomonitor(mi_prin)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#'
#'
#' # averaging
#' data_ts_av <- average_traits(data_ts)
#'
#' numerator <- c("CURRENT_3", "CURRENT_4")
#' denominator <- c("CURRENT_1", "CURRENT_2", "CURRENT_3", "CURRENT_4")
#' fuzzy_trait_ratio(data_agr, tax_lev = "Taxa", trait_db = data_ts_av,
#'    type = "ab", numerator = numerator, denominator = denominator, trans = log1p)


fuzzy_trait_ratio <- function(x, tax_lev = "Taxa", trait_db = NULL, type = "ab", numerator = NULL, denominator = NULL, trans = NULL){

  # useful for transforming data to 0-1 later
  if (inherits(x, "bin")) {
    BIN <- TRUE
  } else {
    BIN <- FALSE
  }

  # set the taxonomic and the trait database
  taxo <- x[[tax_lev]]
  traits <- as.data.frame(trait_db)

  # set the first column of both dataset to Taxa for further merging
  colnames(taxo)[1] <- colnames(traits)[1] <- "Taxa"

  common_names <- intersect(taxo[, "Taxa"], traits[, "Taxa"])

  taxo <- taxo[taxo[, "Taxa"] %in% common_names, ]
  taxo <- taxo[order(taxo[, "Taxa"]), ]

  traits <- traits[traits[, "Taxa"] %in% common_names, ]
  traits <- traits[order(traits[, "Taxa"]), ]

  if(any(taxo[, "Taxa"] != traits[, "Taxa"])){
    stop("Something went wrong, contact the developer")
  }

  # subset the trait dataset
  numerator <- apply(traits[, numerator, drop = FALSE], 1, sum)
  denominator <- apply(traits[, denominator, drop = FALSE], 1, sum)

  i <- unlist(lapply(taxo, is.numeric))
  taxo <- taxo[, i, drop = FALSE]

  if(identical(type, "ab")){
    if(! is.null(trans)){
      taxo <- trans(taxo)
    }
  }

  if(identical(type, "po") | inherits(x, "bin")){
    taxo[taxo > 0] <- 1
  }


  numerator <- apply(sweep(taxo, 1, numerator, "*"), 2, sum)
  denominator <- apply(sweep(taxo, 1, denominator, "*"), 2, sum)

  res <- numerator/denominator
  return(res)

  }
