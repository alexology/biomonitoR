#' @name get_iucn_taxa_tree
#' @title Create biomonitoR reference dataset from GBIF
#'
#' @description Function to create a reference biomonitoR dataset using API from the Global Biodiversity Information Facility (GBIF).
#' @param x Vector containing your taxa list.
#' @param ref_from_tree Create a reference database in the `biomonitoR` format. See [ref_from_tree].
#' @param token See IUCN API token. Detailed information at https://apiv3.iucnredlist.org/api/v3/token
#' @seealso [as_biomonitor] [get_gbif_taxa_tree] [get_nbn_taxa_tree] [get_worms_taxa_tree]
#' @export
#' @examples
#' \dontrun{
#'
#' data(macro_ex)
#' refDB <- get_iucn_taxa_tree(macro_ex$Taxa)
#' refDB_bior <- refDB$taxonomy
#' data_asb <- as_biomonitor(macro_ex, dfref = refDB_bior, overwrite = TRUE)
#' data_agR <- aggregate_taxa(data_asb)
#'
#' }


get_iucn_taxa_tree <- function(x, ref_from_tree = FALSE, token = NULL) {

  # Probably IUCN work only with species names.

  tax <- data.frame() # Create dataframe to store biomonitoR taxa
  notFind <- data.frame() # Create dataframe to store not find taxa
  synonym <- data.frame() # Create dataframe containing synonym names
  # multiple <- data.frame()

  for(j in 1:length(x)){
    print(paste("----",j, "of", length(x), "----", x[j], "----"))
    taxon <- gsub(" ", "%20", x[j])
    url.1 <- paste0("https://apiv3.iucnredlist.org/api/v3/species/", taxon,"/?token=",  token) # Url IUCN
    content.1 <- fromJSON(file = url.1)

    # Url for synonyms
    url.syn <- paste0("https://apiv3.iucnredlist.org/api/v3/species/synonym/", taxon,"/?token=",  token)
    content.syn <- fromJSON(file = url.syn)

    #   Example: Fratercula arctica (European assessment)
    #   tryCatch(
    #   {"Try part: define the expression(s) you want to try"
    #     fromJSON(file = url.1)
    #   },
    #   # Handler when an error occurs:
    #   error = function(cond) {
    #     # Choose a return value when such a type of condition occurs
    #     return(NULL)
    #   }
    # )

    if(!is.null(unlist(content.1$result)) && content.syn$count == 0){

        content.1 <- content.1$result[[1]]

        taxa.bior <-
          data.frame(
            Phylum = upperFirst(only_char(om(content.1$phylum))),
            Class = upperFirst(only_char(om(content.1$class))),
            Subclass = NA,
            Order = upperFirst(only_char(om(content.1$order))),
            Family = upperFirst(only_char(om(content.1$family))),
            Subfamily = NA,
            Tribus = NA,
            Genus = only_char(om(content.1$genus)),
            Species = specieName(only_char(om(content.1$scientific_name))) ,
            Subspecies = if(specieName(only_char(om(content.1$scientific_name))) == 3){
              content.2$valid_name
            } else { NA },
            Taxa = content.1$scientific_name)

        tax <- rbind(tax, taxa.bior)


    } else {
      notFind.1 <- data.frame(taxa = x[j])
      notFind <- rbind(notFind, notFind.1)
    }

    # If there are synonyms
    if(is.null(unlist(content.1$result)) && content.syn$count != 0){

      taxon.2 <- gsub(" ", "%20", content.syn$result[[1]]$accepted_name)
      url.2 <- paste0("https://apiv3.iucnredlist.org/api/v3/species/", taxon.2,
                      "/?token=",  token) # Url IUCN if exists a synonym
      content.2 <- fromJSON(file = url.2)

      content.2 <- content.2$result[[1]]

      taxa.bior <-
        data.frame(
          Phylum = upperFirst(only_char(om(content.2$phylum))),
          Class = upperFirst(only_char(om(content.2$class))),
          Subclass = NA,
          Order = upperFirst(only_char(om(content.2$order))),
          Family = upperFirst(only_char(om(content.2$family))),
          Subfamily = NA,
          Tribus = NA,
          Genus = only_char(om(content.2$genus)),
          Species = specieName(only_char(om(content.2$scientific_name))) ,
          Subspecies = if(specieName(only_char(om(content.2$scientific_name))) == 3){
            content.2$valid_name
          } else { NA },
          Taxa = content.2$scientific_name)

      tax <- rbind(tax, taxa.bior)

      synonym.1 <- data.frame(synonym = x[j],
                              accepted = content.2$scientific_name)
      synonym <- rbind(synonym.1, synonym)
    }


    # if(!is.null(content.1) && lengths(content.1) > 1) {
    #   multiple.1 <- data.frame(taxa = x[j])
    #   multiple <- rbind(multiple, multiple.1)
    # }
  }

  # Message
  if(nrow(notFind) >= 1){
    print(paste("Taxa not found:", nrow(notFind)))
  }
  if(nrow(synonym) >= 1){
    print(paste("Synonym detected:", nrow(synonym)))
  }
  # if(nrow(multiple) >= 1){
  #   print(paste("Taxa with multiple match detected:", nrow(multiple)))
  # }

  if(isTRUE(ref_from_tree)) {
    return(list(
      taxonomy = ref_from_tree(tax[ ,1:ncol(tax)-1]),
      notFindTaxa = notFind,
      synonym = synonym
      #multiMatch = multiple
    ))
  } else {
    return(list(
      taxonomy = tax,
      notFindTaxa = notFind,
      synonym = synonym
      #multiMatch = multiple
    ))

  }
}
