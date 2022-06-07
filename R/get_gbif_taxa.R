#' @name get_gbif_taxa_tree
#' @title Create `biomonitoR` dataset from GBIF
#'
#' @description Function to create a reference dataset using API from the Global Biodiversity Information Facility (GBIF).
#' @param x Vector containing your taxa list.
#' @param ref_from_tree Create a reference database in the `biomonitoR` format. See [ref_from_tree].
#' @keywords as_biomonitor
#' @seealso [as_biomonitor]
#' @importFrom rjson fromJSON
#' @export
#' @examples
#' \dontrun{
#'
#' data(macro_ex)
#' dfref_gbif <- get_gbif_taxa_tree(macro_ex[, "Taxa"])
#' data_asb <- as_biomonitor(macro_ex, dfref = dfref_gbif$taxonomy)
#' data_agg <- aggregate_taxa(data_asb)
#' }



get_gbif_taxa_tree <- function(x, ref_from_tree = FALSE) {

  tax <- data.frame() # Create dataframe to store biomonitoR taxa
  notFind <- data.frame() # Create dataframe to store not find taxa
  synonym <- data.frame() # Create dataframe containing synonym names

  for(j in 1:length(x)){
    print(paste("----",j, "of", length(x), "----", x[j], "----"))
    taxon <- gsub(" ", "%20", x[j])
    url.1 <- paste0("http://api.gbif.org/v1/species/match?&name=", taxon) # GBIF API
    content.1 <- fromJSON(file = url.1)

    if(content.1$matchType != "NONE"){
      if(isTRUE(content.1$synonym)){
        key <- content.1$acceptedUsageKey
        # If co content.1$synonym search the accepted name
        url.2 <- paste0("https://api.gbif.org/v1/species/", key)
        content.2 <- fromJSON(file = url.2)

        synonym.1 <- data.frame(synonym = x[j],
                                accepted = content.2$canonicalName)
        synonym <- rbind(synonym.1, synonym)

      } else {
        key <- content.1$usageKey
        content.2 <- content.1
      }

      taxa.bior <-
        data.frame(
          Phylum = only_char(om(content.2$phylum)),
          Class = only_char(om(content.2$class)),
          Subclass = NA,
          Order = only_char(om(content.2$order)),
          Family = only_char(om(content.2$family)),
          Subfamily = NA,
          Tribus = NA,
          Genus = only_char(om(content.2$genus)),
          Species = only_char(om(content.2$species)), # Is synonym = FALSE I want the NOT accepted name
          Subspecies = ifelse(content.2$rank == "SUBSPECIES",
                              gsub("subsp. ", "", content.1$scientificName),
                              NA),
          Taxa = content.2$canonicalName)

      tax <- rbind(tax, taxa.bior)

    } else {
      notFind.1 <- data.frame(taxa = x[j])
      notFind <- rbind(notFind, notFind.1)
    }
  }

  # Message
  if(nrow(notFind) >= 1){
    print(paste("Taxa not found:", nrow(notFind)))
  }
  if(nrow(synonym) >= 1){
    print(paste("Synonym detected:", nrow(synonym)))
  }

  if(isTRUE(ref_from_tree)) {
    return(list(
      taxonomy = ref_from_tree(tax[ ,1:ncol(tax)-1]),
      notFindTaxa = notFind,
      synonym = synonym
    ))
  } else {
    return(list(
      taxonomy = tax,
      notFindTaxa = notFind,
      synonym = synonym
    ))
  }
}
