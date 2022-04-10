#' @name get_gbif_taxa_tree
#' @title Create biomonitoR dataset from GBIF
#'
#' @description Function to create a reference dataset using API from the Global Biodiversity Information Facility (GBIF).
#' @param x Vector containing your taxa list.
#' @param synonym If `TRUE` synonyms are substituted with their accepted name.
#' @param ref_from_tree Create a reference database in the `biomonitoR` format. See [ref_from_tree].
#' @keywords as_biomonitor
#' @seealso [as_biomonitor]
#' @importFrom rjson fromJSON
#' @export
#' @examples
#' \dontrun{
#'
#' dfref_gbif <- get_gbif_taxa_tree(macro_ex[, "Taxa"])
#' data_asb <- as_biomonitor(macro_ex, dfref = dfref_gbif$taxonomy)
#' data_agg <- aggregate_taxa(data_asb)
#' } data(macro_ex)



get_gbif_taxa_tree <- function(x, synonym = FALSE, ref_from_tree = FALSE) {
  tax <- data.frame() # Create dataframe to store biomonitoR taxa
  statistic <- data.frame() # Create dataframe for statistics
  notFind <- data.frame() # Create dataframe to store not find taxa

  # If SYNONYM == FALSE
  if (isFALSE(synonym)) {
    for (i in 1:length(x)) {
      # print(paste("----",i, "of", length(x), "----", x[i], "----"))
      taxon <- gsub(" ", "%20", x[i])
      ur.1 <- paste0("http://api.gbif.org/v1/species/match?&name=", taxon) # GBIF API
      content.1 <- fromJSON(file = ur.1)

      if (content.1$matchType != "NONE") {
        taxa.bior <-
          data.frame(
            Phylum = only_char(om(content.1$phylum)),
            Class = only_char(om(content.1$class)),
            Subclass = NA,
            Order = only_char(om(content.1$order)),
            Family = only_char(om(content.1$family)),
            Subfamily = NA,
            Tribus = NA,
            Genus = only_char(om(content.1$genus)),
            Species = only_char(om(content.1$species)), # Is synonym = FALSE I want the NOT accepted name
            Subspecies = ifelse(content.1$rank == "SUBSPECIES",
              gsub("subsp. ", "", content.1$scientificName),
              NA
            ),
            # ifelse(content.1$rank == "SUBSPECIES", content.1$scientificName, NA),
            Taxa = content.1$canonicalName
          )

        statistic.bior <-
          data.frame(
            original_taxa = x[i],
            accepted_name = only_char(om(content.1$species)),
            rank = only_char(om(content.1$rank)),
            status = only_char(om(content.1$status)),
            confidence = only_char(om(content.1$confidence))
          )

        tax <- rbind(tax, taxa.bior)
        statistic <- rbind(statistic, statistic.bior)
      } else {
        notFind.1 <- data.frame(taxa = x[i])
        notFind <- rbind(notFind, notFind.1)
      }
    }
  }

  # If SYNONYM == TRUE
  if (isTRUE(synonym)) {
    for (i in 1:length(x)) {
      # print(paste("----",i, "of", length(x), "----", x[i], "----"))
      taxon <- gsub(" ", "%20", x[i])
      ur.1 <- paste0("http://api.gbif.org/v1/species/match?&name=", taxon)
      content.1 <- fromJSON(file = ur.1)

      if (content.1$matchType != "NONE") {
        if (content.1$status != "SYNONYM") {
          taxa.bior <-
            data.frame(
              Phylum = only_char(om(content.1$phylum)),
              Class = only_char(om(content.1$class)),
              Subclass = NA,
              Order = only_char(om(content.1$order)),
              Family = only_char(om(content.1$family)),
              Subfamily = NA,
              Tribus = NA,
              Genus = only_char(om(content.1$genus)),
              Species = only_char(om(content.1$species)),
              Subspecies = ifelse(content.1$rank == "SUBSPECIES",
                gsub("subsp. ", "", content.1$scientificName),
                NA
              ),
              Taxa = content.1$canonicalName
            )

          statistic.bior <-
            data.frame(
              original_taxa = x[i],
              accepted_name = only_char(om(content.1$species)),
              rank = only_char(om(content.1$rank)),
              status = only_char(om(content.1$status)),
              confidence = only_char(om(content.1$confidence))
            )

          tax <- rbind(tax, taxa.bior)
          statistic <- rbind(statistic, statistic.bior)
        }

        if (content.1$status == "SYNONYM") {
          usageKey <- content.1$speciesKey
          ur.2 <- paste0("http://api.gbif.org/v1/species/", usageKey)
          content.2 <- fromJSON(file = ur.2)

          taxa.bior <-
            data.frame(
              Phylum = only_char(om(content.1$phylum)),
              Class = only_char(om(content.1$class)),
              Subclass = NA,
              Order = only_char(om(content.1$order)),
              Family = only_char(om(content.1$family)),
              Subfamily = NA,
              Tribus = NA,
              Genus = only_char(om(content.1$genus)),
              Species = only_char(om(content.1$species)), # If synonym = TRUE I want the accepted name
              Subspecies = ifelse(content.1$rank == "SUBSPECIES",
                gsub("subsp. ", "", content.1$scientificName),
                NA
              ),
              Taxa = content.1$species
            )

          statistic.bior <-
            data.frame(
              original_taxa = x[i],
              accepted_name = only_char(om(content.1$species)),
              rank = only_char(om(content.1$rank)),
              status = only_char(om(content.1$status)),
              confidence = only_char(om(content.1$confidence))
            )

          tax <- rbind(tax, taxa.bior)
          statistic <- rbind(statistic, statistic.bior)
        }
      } else {
        notFind.1 <- data.frame(taxa = x[i])
        notFind <- rbind(notFind, notFind.1)
      }
    }
  }

  tax <- tax[!duplicated(tax$Taxa), ] # Check is there are duplicated names
  tax[is.na(tax)] <- ""

  # # Print the summary results
  # print(list(
  #     rank = table(statistic$rank),
  #     status = table(statistic$status),
  #     notFInd = notFind
  #     )
  #   )

  if (isTRUE(ref_from_tree)) {
    return(list(
      taxonomy = ref_from_tree(tax[, 1:ncol(tax) - 1]),
      taxStat = statistic,
      notFindTaxa = notFind
    ))
  } else {
    return(list(
      taxonomy = tax,
      taxStat = statistic,
      notFindTaxa = notFind
    ))
  }
}
