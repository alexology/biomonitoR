#' Community trait specialization
#'
#' This function calculates the community trait specialization.
#'
#' This function first takes the abundance table corresponding to the desired
#' taxonomic level from the `x` aggregatoR object.
#'
#' Then it searches from the trait data base all the information available at
#' the desired level and, if required, calculates the corresponding averaged
#' trait values (e.g. the family trait values are obtained by averaging all the
#' trait values from taxa with trait information within this family).
#'
#' For each taxon and each trait, a taxon specialization index is calculated
#' using the following formula:
#'
#' `TSI = (sum(c²_tik) - 1/K) / (1 - 1/K)`
#'
#' with `c²_tik` being the affinity of taxon `t` for the modality `k` (among
#' `K`) of trait `i`.
#'
#' Finally, the community trait specialization is calculated for each trait by
#' averaging the TSIs using the transformed abundances (using the `trans`
#' function) as weigths.
#'
#' @inheritParams cwm
#'
#' @return a table with the CSI values for each trait
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#'   n_distinct
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#' data(traitsTachet)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#'
#' csi(x = data.agR, taxLev = "Taxa", trans = log1p)
#' csi(x = data.agR, taxLev = "Taxa",
#'     trans = function(x) {
#'         ifelse(x > 0, 1, 0)
#'     })
#' csi(x = data.agR, taxLev = "Genus", trans = log1p)
#'
#' @seealso [aggregatoR]
#'
#' @references Mondy, C. P., & Usseglio‐Polatera P. (2013) Using Fuzzy-Coded
#'   Traits to Elucidate the Non-Random Role of Anthropogenic Stress in the
#'   Functional Homogenisation of Invertebrate Assemblages. Freshwater Biology,
#'   59 (3), 584‑600. https://doi.org/10.1111/fwb.12289.
#' @references Tachet, H., Richoux, P., Bournaud, M., & Usseglio-Polatera, P.
#'   (2010). Invertébrés d'eau douce: systématique, biologie, écologie. Paris:
#'   CNRS editions.
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www. freshwaterecology.
#'   info–An online tool that unifies, standardises and codifies more than
#'   20,000 European freshwater organisms and their ecological preferences.
#'   Ecological indicators, 53, 271-282.
#'
#' @export

csi <- function(x, traitDB = NULL, taxLev = "Taxa", trans = log1p) {

  # create dummy variables to avoid R CMD check NOTES
  Taxa <- Trait <- Modality <-  Abundance <- Sample <- Weight <- totWeight <-
    k <- TSI <- CSI <- weightedTSI <-  . <- NULL


  tsi <- assignTraits(x = x, traitDB = traitDB, taxLev = taxLev) %>%
    gather(key = Modality, value = Affinity, -Taxa)              %>%
    mutate(Trait = strsplit(Modality, split = "_")               %>%
             sapply(FUN = '[[', 1))                              %>%
    group_by(Taxa, Trait)                                        %>%
    mutate(k = n_distinct(Modality))                             %>%
    summarise(TSI = (sum(Affinity^2) - 1 / unique(k)) / (1 - 1 / unique(k)))

  abundances <- x[[taxLev]]
  colnames(abundances)[1] <- "Taxa"

  abundances                                         %>%
    gather(key = Sample, value = Abundance, -Taxa)   %>%
    mutate(Sample = factor(Sample,
                           levels = colnames(abundances)[-1]),
           Taxa   = as.character(Taxa),
           Weight = trans(Abundance))                %>%
    left_join(group_by(., Sample)                    %>%
                summarise(totWeight = sum(Weight)),
              by = "Sample")                         %>%
    left_join(tsi,
              by = "Taxa")                           %>%
    mutate(weightedTSI = (Weight * TSI) / totWeight) %>%
    group_by(Sample, Trait)                          %>%
    summarise(CSI = sum(weightedTSI, na.rm = TRUE))  %>%
    spread(key = Trait, value = CSI)                 %>%
    as.data.frame()

}
