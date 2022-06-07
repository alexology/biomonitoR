#' @name get_worms_taxa_tree
#' @title Function for creating reference databases using API from World Register of Marine Species (WORMS)
#'
#' @description This function allows creating our custom reference database using API from World Register of Marine Species (WORMS).
#' @param x Vector contain your Taxa list.
#' @param ref_from_tree Create a reference database in the `biomonitoR` format. See [ref_from_tree].
#' @seealso [as_biomonitor]
#' @export
#' @examples
#' \dontrun{
#'
#' data(macro_ex)
#' dfref_worms <- get_worms_taxa_tree(macro_ex[, "Taxa"])
#' data_asb <- as_biomonitor(macro_ex, dfref = dfref_worms$taxonomy)
#' data_agg <- aggregate_taxa(data_asb)
#'
#' }

get_worms_taxa_tree <- function(x, ref_from_tree = FALSE) {
  tax <- data.frame() # Create dataframe to store biomonitoR taxa
  notFind <- data.frame() # Create dataframe to store not find taxa
  synonym <- data.frame() # Create dataframe containing synonym names
  multiple <- data.frame() # Create dataframe containing taxa with multiple matches (-999)

  for (j in 1:length(x)) {
    print(paste("----", j, "of", length(x), "----", x[j], "----"))
    taxon <- gsub(" ", "%20", x[j])
    url.aphia <- paste0("https://www.marinespecies.org/rest/AphiaIDByName/", taxon, "?marine_only=false")
    aphia <- tryCatch(
      {
        "Try part: define the expression(s) you want to try"
        fromJSON(file = url.aphia) # Retrieve Aphia
      },
      # Handler when an error occurs:
      error = function(cond) {
        # Choose a return value when such a type of condition occurs
        return(NULL)
      }
    )

    if (!is.null(aphia) && aphia != -999) {
      # Retrieve basic information about the taxon
      url.1 <- paste0("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/", aphia)
      content.1 <- fromJSON(file = url.1)
      if (content.1$status != "accepted") {
        key <- content.1$valid_AphiaID
        # If content.1$status != "accepted" search the accepted name
        url.2 <- paste0("https://www.marinespecies.org/rest/AphiaRecordByAphiaID/", key)
        content.2 <- fromJSON(file = url.2)

        synonym.1 <- data.frame(
          synonym = x[j],
          accepted = content.2$valid_name
        )
        synonym <- rbind(synonym.1, synonym)
      } else {
        # key <- content.1$AphiaID
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
          Species = specieWorms(only_char(om(content.2$valid_name))),
          Subspecies = if (countWords(only_char(om(content.2$valid_name))) == 3) {
            content.2$valid_name
          } else {
            NA
          },
          Taxa = content.2$valid_name
        )

      tax <- rbind(tax, taxa.bior)
    }

    if (is.null(aphia)) {
      notFind.1 <- data.frame(taxa = x[j])
      notFind <- rbind(notFind, notFind.1)
    }

    if (!is.null(aphia) && aphia == -999) {
      multiple.1 <- data.frame(taxa = x[j])
      multiple <- rbind(multiple, multiple.1)
    }
  }

  # Message
  if (nrow(notFind) >= 1) {
    print(paste("Taxa not found:", nrow(notFind)))
  }
  if (nrow(synonym) >= 1) {
    print(paste("Synonym detected:", nrow(synonym)))
  }
  if (nrow(multiple) >= 1) {
    print(paste("Taxa with multiple match detected:", nrow(multiple)))
  }


  if (isTRUE(ref_from_tree)) {
    return(list(
      taxonomy = ref_from_tree(tax[, 1:ncol(tax) - 1]),
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
