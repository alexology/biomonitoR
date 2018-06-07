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
#'   By default, the data base used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
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
#' data(traitsTachet)
#'
#' data.bio <- asBiomonitor(macro_ex)
#' data.agR <- aggregatoR(data.bio)
#'
#' cwm(x = data.agR, traitDB = traitsTachet, taxLev = "Taxa", trans = log1p)
#' cwm(x = data.agR, traitDB = traitsTachet, taxLev = "Taxa",
#'     trans = function(x) {
#'         ifelse(x > 0, 1, 0)
#'     })
#' cwm(x = data.agR, traitDB = traitsTachet, taxLev = "Genus", trans = log1p)
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

cwm <- function(x, traitDB = "traitsTachet", taxLev = "Taxa", trans = log1p) {

  # create dummy variables to avoid R CMD check NOTES
  traitsTachet <- Taxa <- modality <- affinity <- Phylum <- Subspecies <-
    Abundance <- Sample <- Weight <- Affinity <- totWeight <-
    weightedAffinity <- Category <- . <- NULL

  # check if the object x is of class "biomonitoR"
  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  if (! taxLev %in% c("Family", "Genus", "Species", "Taxa")) {
    return("taxLev should be one of the following: Family, Genus, Species or Taxa")
  }


  trait_db <- traitDB                               %>%
    (function(df) {
      mutate(df,
             Taxa = gsub(pattern     = " sp.",
                          replacement = "",
                          x           = Taxa))
    })                                              %>%
    gather(key = modality, value = affinity, -Taxa) %>%
    group_by(Taxa, modality)                        %>%
    summarise(affinity = mean(affinity))            %>%
    spread(key = modality, value = affinity)        %>%
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

  taxa_traits <- mutate(ref, Taxa = as.character(Taxa)) %>%
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
    })


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
