#' @name get_nbn_taxa_tree
#' @title Create `biomonitoR` reference dataset from the National Biodiversity Network (NBN)
#'
#' @description Function to create a reference `biomonitoR` dataset using API from the National Biodiversity Network (NBN).
#' @param x Vector contain your Taxa list. This is the base list to create your custom reference database.
#' @param ref_from_tree Create a reference database in the `biomonitoR` format. See [ref_from_tree].
#' @seealso [as_biomonitor] [get_gbif_taxa_tree] [get_iucn_taxa_tree] [get_worms_taxa_tree]
#' @export
#' @examples
#'
#' \dontrun{
#'
#' data(macro_ex)
#' dfref_c <- get_gbif_taxa_tree(macro_ex[, "Taxa"])
#' data_asb <- as_biomonitor(macro_ex, dfref = dfref_nbn$taxonomy)
#' data_agg <- aggregate_taxa(data_asb)
#'
#' }



get_nbn_taxa_tree <- function(x, ref_from_tree = FALSE) {

  tax <- data.frame() # Create dataframe to store biomonitoR taxa
  notFind <- data.frame() # Create dataframe to store not find taxa
  synonym <- data.frame() # Create dataframe containing synonym names
  multiple <- data.frame() # Create dataframe containing taxa with multiple matches (-999)

  #URL 1. Search ID
  # https://species-ws.nbnatlas.org/guid/batch?q=Fox

  # URL 2. Taxonomy
  # https://species-ws.nbnatlas.org/classification/NBNSYS0000005616

  # In this case to verify the synonyms we can check the accordance between the fields identifier and acceptedIdentifier. If they are equal, the name is NOT a synonym.

  for(j in 1:length(x)){
    print(paste("----",j, "of", length(x), "----", x[j], "----"))
    taxon <- gsub(" ", "%20", x[j])
    url.1 <- paste0("https://species-ws.nbnatlas.org/guid/batch?q=", taxon) # Url NBN
    content.1 <- tryCatch(
      {"Try part: define the expression(s) you want to try"
        fromJSON(file = url.1)
      },
      # Handler when an error occurs:
      error = function(cond) {
        # Choose a return value when such a type of condition occurs
        return(NULL)
      }
    )

    if(!is.null(content.1) && lengths(content.1) == 1){
      content.1 <- content.1[[1]][[1]]
      if(isFALSE(content.1$identifier == content.1$acceptedIdentifier)){
        # For synonyms
        key <- content.1$acceptedIdentifier
        # Search accepted name
        url.2 <- paste0("https://species-ws.nbnatlas.org/classification/", key)
        content.2 <- fromJSON(file = url.2)

        synonym.1 <- data.frame(synonym = x[j],
                                accepted = content.2[[10]]$scientificName)
        synonym <- rbind(synonym.1, synonym)

      } else {
        key <- content.1$identifier

        url.2 <- paste0("https://species-ws.nbnatlas.org/classification/", key)
        content.2 <- fromJSON(file = url.2)

      }

      content.2.df <- data.frame(Reduce(rbind, content.2))

      taxa.bior <-
        data.frame(
          Phylum = ifelse("phylum" %in% content.2.df$rank,
                          only_char(om(content.2.df$scientificName[content.2.df$rank == "phylum"][[1]])),
                          NA),
          Class = ifelse("class" %in% content.2.df$rank,
                         only_char(om(content.2.df$scientificName[content.2.df$rank == "class"][[1]])),
                         NA),
          Subclass = NA,
          Order = ifelse("order" %in% content.2.df$rank,
                         only_char(om(content.2.df$scientificName[content.2.df$rank == "order"][[1]])),
                         NA),
          Family = ifelse("family" %in% content.2.df$rank,
                          only_char(om(content.2.df$scientificName[content.2.df$rank == "family"][[1]])),
                          NA),
          Subfamily = NA,
          Tribus = NA,
          Genus = ifelse("genus" %in% content.2.df$rank,
                         only_char(om(content.2.df$scientificName[content.2.df$rank == "genus"][[1]])),
                         NA),
          Species = ifelse("species" %in% content.2.df$rank,
                           only_char(om(content.2.df$scientificName[content.2.df$rank == "species"][[1]])),
                           NA),
          Subspecies = ifelse("subspecies" %in% content.2.df$rank,
                              only_char(om(content.2.df$scientificName[content.2.df$rank == "subspecies"][[1]])),
                              NA),
          Taxa = as.character(content.2.df$scientificName[nrow(content.2.df)]))

      tax <- rbind(tax, taxa.bior)

    }

    # Taxa with multiple matches
    if(!is.null(content.1) && lengths(content.1) > 1) {
      multiple.1 <- data.frame(taxa = x[j])
      multiple <- rbind(multiple, multiple.1)

      content.1.df <- t(as.data.frame(content.1[[1]]))
      content.1.name <- content.1.df[grepl("acceptedName", rownames(content.1.df)), ]
      name.selection <- menu(c(content.1.name), title="Chose options")

      # Accepted key of the selected taxon
      key <- content.1[[1]][[name.selection]]$acceptedIdentifier
      # Search selected name from multiplenames
      url.3 <- paste0("https://species-ws.nbnatlas.org/classification/", key)
      content.3 <- fromJSON(file = url.3)
      content.3.df <- data.frame(Reduce(rbind, content.3))

      taxa.bior <-
        data.frame(
          Phylum = ifelse("phylum" %in% content.3.df$rank,
                          only_char(om(content.3.df$scientificName[content.3.df$rank == "phylum"][[1]])),
                          NA),
          Class = ifelse("class" %in% content.3.df$rank,
                         only_char(om(content.3.df$scientificName[content.3.df$rank == "class"][[1]])),
                         NA),
          Subclass = NA,
          Order = ifelse("order" %in% content.3.df$rank,
                         only_char(om(content.3.df$scientificName[content.3.df$rank == "order"][[1]])),
                         NA),
          Family = ifelse("family" %in% content.3.df$rank,
                          only_char(om(content.3.df$scientificName[content.3.df$rank == "family"][[1]])),
                          NA),
          Subfamily = NA,
          Tribus = NA,
          Genus = ifelse("genus" %in% content.3.df$rank,
                         only_char(om(content.3.df$scientificName[content.3.df$rank == "genus"][[1]])),
                         NA),
          Species = ifelse("species" %in% content.3.df$rank,
                           only_char(om(content.3.df$scientificName[content.3.df$rank == "species"][[1]])),
                           NA),
          Subspecies = ifelse("subspecies" %in% content.3.df$rank,
                              only_char(om(content.3.df$scientificName[content.3.df$rank == "subspecies"][[1]])),
                              NA),
          Taxa = as.character(content.3.df$scientificName[nrow(content.3.df)])
        )

      tax <- rbind(tax, taxa.bior)

    }

    if(isFALSE(!is.null(content.1) && lengths(content.1) == 1) &&
       isFALSE(!is.null(content.1) && lengths(content.1) > 1)){
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
  if(nrow(multiple) >= 1){
    print(paste("Taxa with multiple match detected:", nrow(multiple)))
  }

  if(isTRUE(ref_from_tree)) {
    return(list(
      taxonomy = ref_from_tree(tax[ ,1:ncol(tax)-1]),
      notFindTaxa = notFind,
      synonym = synonym,
      multiMatch = multiple
    ))
  } else {
    return(list(
      taxonomy = tax,
      notFindTaxa = notFind,
      synonym = synonym,
      multiMatch = multiple
    ))

  }
}
