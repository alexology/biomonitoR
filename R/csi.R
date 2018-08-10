#' Community trait specialization
#'
#' This function calculates the community trait specialization.
#' @param x results of function aggregatoR
#' @param traitDB a trait data base with a column `Taxa` and the other columns
#'   containing the traits. If a trait has several modalities they should be
#'   named as follow: TRAIT_MODALITY.
#'
#'   By default, the data base used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#' @param taxLev character string giving the taxonomic level used to retrieve
#'   trait information. Possible levels are `"Taxa"`, `"Species"`, `"Genus"`,
#'   `"Family"` as returned by the [aggregatoR] function.
#' @param trans the function used to transform the abundances, by default
#'   [log1p]
#' @details This function first takes the abundance table corresponding to the desired
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
#' \deqn{TSI = (sum(c^{2}_tik) - 1/K) / (1 - 1/K)}
#'
#' with \eqn{c^{2}_tik} being the affinity of taxon `t` for the modality `k` (among
#' `K`) of trait `i`.
#'
#' Finally, the community trait specialization is calculated for each trait by
#' averaging the TSIs using the transformed abundances (using the `trans`
#' function) as weigths.
#'
#'
#' @return a table with the CSI values for each trait
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup n_distinct
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#' data(traitsTachet)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#'
#' csi(x = data.agR, traitDB = NULL, taxLev = "Taxa", trans = log1p)
#' csi(x = data.agR, traitDB = NULL, taxLev = "Taxa",
#'     trans = function(x) {
#'         ifelse(x > 0, 1, 0)
#'     })
#' csi(x = data.agR, traitDB = NULL, taxLev = "Genus", trans = log1p)
#'
#' @seealso [aggregatoR]
#'
#' @references Mondy, C. P., & Usseglio-Polatera P. (2013) Using Fuzzy-Coded
#'   Traits to Elucidate the Non-Random Role of Anthropogenic Stress in the
#'   Functional Homogenisation of Invertebrate Assemblages. Freshwater Biology,
#'   59 (3), 584-600. https://doi.org/10.1111/fwb.12289.
#' @references Tachet, H., Richoux, P., Bournaud, M., & Usseglio-Polatera, P.
#'   (2010). Invertebres d'eau douce: systematique, biologie, ecologie. Paris:
#'   CNRS editions.
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www. freshwaterecology.
#'   info-An online tool that unifies, standardises and codifies more than
#'   20,000 European freshwater organisms and their ecological preferences.
#'   Ecological indicators, 53, 271-282.
#'
#' @export

csi <- function(x, traitDB = NULL, taxLev = "Taxa", trans = log1p) {

  if( is.null( traitDB )){
    traitDB = traitsTachet
  } else {
    traitDB = traitDB
  }

  # create dummy variables to avoid R CMD check NOTES
  Taxa <- Trait <- Modality <- Affinity <- Phylum <- Subspecies <-
    Abundance <- Sample <- Weight <- totWeight <- k <-
    TSI <- CSI <- weightedTSI <-  . <- NULL

  # check if the object x is of class "biomonitoR"
  classCheck(x, group = "mi")

  if (! taxLev %in% c("Family", "Genus", "Species", "Taxa")) {
    return("taxLev should be one of the following: Family, Genus, Species or Taxa")
  }

  trait_db <- traitDB                               %>%
    (function(df) {
      mutate(df,
             Taxa = gsub(pattern     = "sp.|Ad.|Lv.|Gen.",
                         replacement = "",
                         x           = Taxa))
    })                                              %>%
    gather(key = Modality, value = Affinity, -Taxa) %>%
    group_by(Taxa, Modality)                        %>%
    summarise(Affinity = mean(Affinity))            %>%
    spread(key = Modality, value = Affinity)        %>%
    ungroup()

  abundances <- x[[taxLev]]
  colnames(abundances)[1] <- "Taxa"

  taxa       <- as.character(abundances$Taxa)

  if (length(taxa[taxa != "unassigned"]) == 0) {
    return("At least one taxa should be identified at a level compatible with the indicated taxLev")
  }

  if (taxLev == "Taxa") {
    level <- sapply(select(x$Tree, Phylum:Subspecies),
                    function(i) {
                      as.character(i) == as.character(x$Tree$Taxa)
                    }) %>%
      (function(df) {
        colnames(df)[apply(df, MARGIN = 1, which)]
      })
  } else {
    level <- rep(taxLev, length(taxa))
  }

  ref <- select(x$Tree, Phylum:Taxa)

  tsi <- mutate(ref, Taxa = as.character(Taxa))         %>%
    left_join(mutate(trait_db, Taxa = as.character(Taxa)),
              by = "Taxa")                              %>%
    (function(df) {
      lapply(seq(length(taxa)),
             function(i) {
               df[df[, level[i]] == taxa[i],] %>%
                 select(-(Phylum:Taxa))       %>%
                 colMeans(na.rm = TRUE)
             })               %>%
        do.call(what = rbind) %>%
        data.frame(Taxa = taxa, ., stringsAsFactors = FALSE)
    })                                                  %>%
    gather(key = Modality, value = Affinity, -Taxa)     %>%
    mutate(Trait = strsplit(Modality, split = "_")      %>%
             sapply(FUN = '[[', 1))                     %>%
    group_by(Taxa, Trait)                               %>%
    mutate(k = n_distinct(Modality))                    %>%
    summarise(TSI = (sum(Affinity^2) - 1 / unique(k)) /
                (1 - 1 / unique(k)))

  abundances                                        %>%
    gather(key = Sample, value = Abundance, -Taxa)  %>%
    mutate(Sample = factor(Sample,
                           levels = colnames(abundances)[-1]),
           Taxa   = as.character(Taxa),
           Weight = trans(Abundance))               %>%
    left_join(group_by(., Sample)                   %>%
                summarise(totWeight = sum(Weight)),
              by = "Sample")                        %>%
    left_join(tsi,
              by = "Taxa")                          %>%
    mutate(weightedTSI = (Weight * TSI) /
             totWeight)                             %>%
    group_by(Sample, Trait)                         %>%
    summarise(CSI = sum(weightedTSI, na.rm = TRUE)) %>%
    spread(key = Trait, value = CSI)                %>%
    as.data.frame()

}
