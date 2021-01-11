#' @describeIn csi taxon specialization index
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup n_distinct semi_join
#' @importFrom tidyr gather spread

tsi <- function(x, trait_db = NULL, tax_lev = "Taxa") {
  classCheck(x)

  if (is.null(trait_db)) {
    trait_db <- traitsTachet
  } else {
    trait_db <- trait_db
  }

  # create dummy variables to avoid R CMD check NOTES

  Taxa <- Trait <- Modality <- modality <- affinity <- Affinity <- Phylum <- Subspecies <-
    Abundance <- Sample <- Weight <- totWeight <- k <-
    TSI <- CSI <- weightedTSI <- . <- NULL


  abundances <- x[[tax_lev]]
  colnames(abundances)[1] <- "Taxa"
  abundances$Taxa <- as.character(abundances$Taxa)

  if (inherits(x, "bin")) {
    abundances <- to_bin(abundances)
  }

  abundances <- merge(abundances, trait_db[, "Taxa", drop = FALSE])
  trait_db <- merge(trait_db, abundances[, "Taxa", drop = FALSE])

  tsi <- trait_db %>%
    mutate(Taxa = as.character(Taxa)) %>%
    semi_join(abundances, "Taxa") %>%
    gather(key = Modality, value = Affinity, -Taxa) %>%
    mutate(Trait = strsplit(Modality, split = "_") %>%
      sapply(FUN = "[[", 1)) %>%
    group_by(Taxa, Trait) %>%
    mutate(k = n_distinct(Modality)) %>%
    summarise(TSI = (sum(Affinity^2) - 1 / unique(k)) /
      (1 - 1 / unique(k))) %>%
    spread(key = Trait, value = TSI) %>%
    as.data.frame()

  tsi
}
