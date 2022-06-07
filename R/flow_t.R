#' @name flow_t
#' @title Flow-T index
#'
#' @description This function calculates the Flow-T index.
#' @param x Result of `aggregate_taxa()`.
#' @param tax_lev Taxonomic level at which the calculation has to be performed.
#' @param type Presence only `po` or abundance `ab`.
#' @param trans The function used to transform the abundances, by default `log1p`.
#'
#' @details This function calculates the Flow_T index developed by Laini et al. (2022).
#' Flow-T is calculated using the current velocity preference categories for each sample using:
#'
#' \deqn{Flow-T= \frac{\sum_{i=1}^{n}\log(A_i + 1)(m_i + f_i)}{\sum_{i=1}^{n}\log(A_i + 1)}}
#'
#' Where n represents the number of taxa in a sample, A_i the abundance of the ith taxon, mi and fi the
#' "moderate" and "fast" velocity preference classes of the ith taxon according to Tachet et al. (2010).
#' The logarithm of the abundance can be replaced by 1 to obtain the presence-absence version of the Flow-T
#' index.
#'
#' @seealso [as_biomonitor] [assign_traits]
#' @export
#' @references Laini, A., Burgazzi, G., Chadd, R., England, J., Tziortzis, I., Ventrucci, M., Vezza, P., Viaroli, P., Wood, P.J. & Guareschi, S. (2022). Using invertebrate functional traits to improve flow variability assessment within European rivers. Science of The Total Environment, 832, 155047.
#' @references Tachet, H., Richoux, P., Bournaud, M., Usseglio-Polatera, P., (2010). Invertebres d'eau douce: Systematique, biologie, ecologie, edition revue et augmentee. CNRS EDITIONS, Paris.
#' @examples
#'
#' data(mi_prin)
#' data_bio <- as_biomonitor(mi_prin)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#'



flow_t <- function(x, tax_lev = "Taxa", type = "ab", trans = NULL){

  data_ts <- assign_traits(x)
  data_ts_av <- average_traits(data_ts)

  numerator <- c("CURRENT_3", "CURRENT_4")
  denominator <- c("CURRENT_1", "CURRENT_2", "CURRENT_3", "CURRENT_4")

  res <- fuzzy_trait_ratio(x, tax_lev = tax_lev, trait_db = data_ts_av,
                    type = type, numerator = numerator, denominator = denominator, trans = trans)

  res
}
