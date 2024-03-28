#' Drought Effect of Habitat Loss on Invertebrates
#'
#' @description
#' This function calculates the *Average Score Per Taxon* index following Armitage et al. (1983), Davy-Bowker et al. (2007) and Alba-Tercedor & Sanchez-Ortega (1988) implementations.
#'
#' @param x Result of `aggregate_taxa()`.
#' @param method The implementation of DEHLI needed to calculate ASPT. The only choice is `chadd_2017` (Chadd et al . 2017).
#'  Users can provide their own data.frame (see examples) with a column called *Taxon* and the column of scores called *Scores*.
#' @param agg This option allows the composite family approach. It can be `FALSE`, `TRUE` or a `data.frame`.
#' If `FALSE` no aggregation will be performed, while if `TRUE` aggregation will be performed according to the rules described in *Details*.
#' A `data.frame` containing the aggregation rules can be provided by the user.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' `agg` cannot be `TRUE` when a `data.frame` is provided as method.
#' @param exceptions Taxa that need to be excluded from the calculation.
#' This option can be useful, for instance, to exclude an alien species belonging to an autochthonous family.
#' @param traceB If set to `TRUE` a list as specified below will be returned.
#'
#'
#' @keywords dehli
#'
#' @references Chadd, R. P., England, J. A., Constable, D., Dunbar, M. J., Extence, C. A., Leeming, D. J., Murray-Bligh Wood J.A., P. J. (2017). An index to track the ecological effects of drought development and recovery on riverine invertebrate communities. Ecological Indicators, 82, 344-356.
#'
#' @details The DEHLI index measures the effect of flow intermittency on the macroinvertebrate community.
#' It assigning weights to taxa on the basis of their likely association with key stages of channel drying.
#' DEHLI is the average score of the taxa receiving a score, similary to [aspt].
#'
#' The follwing composite taxa are used:
#' \enumerate{
#'   \item Philopotamus (inc. Wormaldia)
#'   \item Brachyptera (inc. Rhabdiopteryx)
#'   \item Protonemura (inc. Amphinemura)
#'   \item Protonemura (inc. Nemurella)
#'   \item Habrophlebia (inc. Leptophlebia)
#' }
#'
#'
#' The `dehli()` function automatically check for parent-child pairs in the scoring system, see the Value section for a definition.
#' All the information used for `dehli()` calculation can be retrieved with the function \code{\link{show_scores}}.
#'
#' @return If `traceB` is set to `TRUE` a list with the following elements will be returned:
#' \itemize{
#'  \item `results` Result of `dehli()`.
#'  \item `taxa_df` The `data.frame` used for the calculation containing the abundance of taxa receiving a score.
#'  \item `composite_taxa` Taxa aggregated following the aggregation rules when `agg` is not `NULL`.
#'  \item `exceptions` A `data.frame` containing the changes made by excluding the taxa included in `exceptions`.
#'  \item `parent_child_pairs` For instance in Spanish ASPT both *Ferrissia* and Planorbidae receive a score.
#'  Abundances of the higher taxonomic level need therefore to be adjusted by subtracting the abundances of the lower taxonomic level.
#' }
#'
#' @seealso [aggregate_taxa] [bmwp]
#' @export
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' dehli(data_agr)
#' dehli(data_agr)


dehli <- function(x, method = "chadd_2017", agg = TRUE, exceptions = NULL, traceB = FALSE){

  if (! identical(method, "chadd_2017") & !(is.data.frame(method))) (stop("Please provide a valid method"))

  if(identical(method, "chadd_2017")){
    method <- dehli_scores_chadd_2017
  } else {
    method <- method
  }

  if(isTRUE(agg)){
    agg <- dehli_acc_chadd_2017
  }

  aspt(x = x, method = method, agg = agg, exceptions = exceptions, traceB = traceB)
}
