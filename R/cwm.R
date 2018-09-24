#' Community-Weighted Mean values of traits
#'
#' This function calculates the community-weighted means of trait categories.
#'
#' This function first takes the abundance table corresponding to the desired
#' taxonomic level from the `x` aggregatoR object.
#'
#' Then it searches from the trait data base all the information available at
#' the desired level and, if required, calculates the corresponding averaged
#' trait values (e.g. the family trait values are obtained by averaging all the
#' trait values from taxa with trait information within this family).
#'
#' Finally, the community mean trait values are calculated using the transformed
#' abundances (using the `trans` function) as weigths
#'
#' @param x results of function aggregatoR
#' @param traitDB a trait data base with a column `Taxa` and the other columns
#'   containing the traits. If a trait has several modalities they should be
#'   named as follow: TRAIT_MODALITY.
#'
#'   By default (NULL), the data base used is the one from Tachet *et al* (2010)
#'   that can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#' @param taxLev character string giving the taxonomic level used to retrieve
#'   trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#'   `"Family"` as returned by the [aggregatoR] function.
#' @param trans the function used to transform the abundances, by default
#'   [log1p]
#'
#' @return a table with the CWM values of each trait (trait modality)
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#'
#' cwm(x = data.agR, taxLev = "Taxa", trans = log1p)
#' cwm(x = data.agR, taxLev = "Taxa",
#'     trans = function(x) {
#'         ifelse(x > 0, 1, 0)
#'     })
#' cwm(x = data.agR, taxLev = "Genus", trans = log1p)
#'
#' @seealso [aggregatoR]
#'
#' @references Tachet, H., Richoux, P., Bournaud, M., & Usseglio-Polatera, P.
#'   (2010). Invertébrés d'eau douce: systématique, biologie, écologie. Paris:
#'   CNRS editions.
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www. freshwaterecology.
#'   info–An online tool that unifies, standardises and codifies more than
#'   20,000 European freshwater organisms and their ecological preferences.
#'   Ecological indicators, 53, 271-282.
#'
#' @export

cwm <- function(x, traitDB = NULL, taxLev = "Taxa", trans = log1p) {

  # create dummy variables to avoid R CMD check NOTES
  Taxa <- Abundance <- Sample <- Weight <- Affinity <- totWeight <-
    weightedAffinity <- Category <- . <- NULL

  taxa_traits <- assignTraits(x = x, traitDB = traitDB, taxLev = taxLev)

  abundances <- x[[taxLev]]
  colnames(abundances)[1] <- "Taxa"

  abundances                                       %>%
    gather(key = Sample, value = Abundance, -Taxa) %>%
    mutate(Sample = factor(Sample,
                           levels = colnames(abundances)[-1]),
           Taxa   = as.character(Taxa),
           Weight = trans(Abundance))              %>%
    left_join(group_by(., Sample)                  %>%
                summarise(totWeight = sum(Weight)),
              by = "Sample")                       %>%
    left_join(gather(taxa_traits,
                     key   = Category,
                     value = Affinity, -Taxa),
              by = "Taxa")                         %>%
    mutate(weightedAffinity = (Weight * Affinity) /
             totWeight)                            %>%
    group_by(Sample, Category)                     %>%
    summarise(Affinity = sum(weightedAffinity,
                             na.rm = TRUE))        %>%
    spread(key = Category, value = Affinity)       %>%
    as.data.frame()

}
