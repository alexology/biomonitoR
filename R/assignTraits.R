#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread

assignTraits <- function(x, traitDB = NULL, taxLev = "Taxa") {
  if( is.null( traitDB )){
    traitDB = traitsTachet
  } else {
    traitDB = traitDB
  }

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
             Taxa = gsub(pattern     = "sp.|Ad.|Lv.|Gen.",
                         replacement = "",
                         x           = Taxa) %>%
               trimws(x = ., which = "both"))
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

  ref <- select(x$Tree, Phylum:Taxa)

  mutate(ref, Taxa = as.character(Taxa)) %>%
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
}
