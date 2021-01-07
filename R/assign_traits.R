#' assign_traits
#'
#' @description
#' A function for scaling traits across taxonomic levels.
#'
#' @details
#' This function allows to obtain missing traits for a target taxon by taking traits from lower or to upper taxomic levels.
#' For instance, consider the case where the genus Acroloxus is present in the user dataset and the species Acroloxus lacustris
#' in the traits database. A simple merge would exclude Acroloxus from the tha analysis since Acroloxus and A. lacustris
#' would not match. The function traitScaling allows to assign Acroloxus lacustris traits to Acroloxus.
#' This function works also in the opposite direction. Consider the case where there are no traits for the target taxon and
#' the target taxon has not been identified at species level. The function traitScaling will assign the traits of the nearest taxonomic level
#' to the target taxa (e.g. Tanypodinae traits assign to Ablabesmyia monilis). Consider also these examples to understand the behaviour of this
#' function. For instance Anabolia lombarda is present in the user taxomic dataset while only Anabolia nervosa and Anabolia are present
#' in the trait database. In this case traitScaling will assing to A. lombarda only the traits of Anabolia.
#' Moreover, let assume that Coelostoma is present in the user dataset while only Berosus and Crenitis punctatostriata are present in the traits database.
#' Here traitScaling will assign to Coelostoma the scores of Berosus and C. punctatostriata since they belong to the same family and there are no information at family level. \cr \cr
#' The function traitScaling will measure the taxonomic distance between the target taxa and the taxa used to assign the trait score. This distance
#' can be negative (e.g. Species to Genus) and positive (e.g. Genus to Species). The distance is measured assigning values as follows:
#' 1 (Species to Genus) , 2 (Species to family), -1 (Family to Genus), etc. traitScaling considers only the
#' taxonomic levels from Subspecies to Family (Subspecies, Species, Genus, Tribus, Subfamily, Family).
#'
#'
#' @param x Results of function `aggregate_taxa()`
#' @param trait_db A trait data base with a column `Taxa` and the other columns
#'   containing the traits.
#'   By default, the database used is the one from Tachet *et al* (2010) that
#'   can be retrieved from
#'   [freshwaterecology.info](https://www.freshwaterecology.info/) website
#'   (Schmidt-Kloiber & Hering, 2015).
#'   It includes traits only for macroinvertebrates.
#' @param group Biotic group of interest. Possible values are `mi` for macroinvertebrates, `mf` for macrophytes and `fi` for fish.
#'  The choice will set the right reference database for the specified group.
#'  This option will not be considered if a custom reference database is provided. Default to `mi`.
#' @param tax_lev Taxonomic level on which the calculation has to be made.
#' Default to `Taxa`, the maximum taxonomic level is `Family`.
#' @param dfref Reference database as used in the function aggregatoR.
#' @param filter_by_distance Filter the results according to the taxonomic distance. Possible values are `pos`, `neg` or a positive integer. See details.
#' @param col_blocks A vector that contains the number of modalities for each trait
#'
#'
#' @importFrom dplyr '%>%' mutate select left_join group_by summarise ungroup
#' @importFrom tidyr gather spread
#'
#' @examples
#' data(macro_ex)
#'
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#'
#' # select only the nearest traits
#' data_ts_sub <- manageTraits(data_ts, method = "nearest+-")
#'
#' # averaging
#' data_ts_av <- average_traits(data_ts_sub)
#'
#' # traits random sampling
#' data_ts_st <- sample_traits(data_ts)
#'
#' @seealso [aggregate_taxa]
#'
#' @export
#'
#' @export sample_traits
#'
#' @export average_traits



assign_traits <- function(x, trait_db = NULL, group = "mi", tax_lev = "Taxa", dfref = NULL, filter_by_distance = NULL) {
  if (is.null(trait_db)) {
    # check if x is of class biomonitoR and mi
    classCheck(x)

    if (inherits(x, "custom") & is.null(trait_db)) {
      warning("It seems that you used your own reference database. Please check the consistency of the taxonomy used for calculating the index with those of your reference database to have reliable results.")
    }

    trait_db <- traitsTachet
  } else {
    trait_db$Taxa <- trimws(trait_db$Taxa)
    classCheck(x)
  }

  if (is.null(dfref)) {
    if (identical(group, "mi")) {
      dfref <- mi_ref
    }
    if (identical(group, "mf")) {
      dfref <- mf_ref
    }
    if (identical(group, "fi")) {
      dfref <- fi_ref
    }
  } else {
    dfref <- dfref
  }



  # merge the trait database with the reference database in order to scale
  # traits across taxonomic levels
  # [ , - 1] deletes the Taxa columns
  ref <- merge(dfref, trait_db, by = "Taxa", sort = FALSE)[, -1]

  # create a data.frame with the same column as trait_db but with 0 rows
  # it will be important later to iterate rbind to this object
  trait.interm <- data.frame(trait_db[-c(1:nrow(trait_db)), ])

  # ref.na is a 0 length vector to store the taxa names of the selected rows
  # this allow to does not loose the information about the name of the taxa at the orginial
  # taxonomic level

  ref.na <- c()

  # 10 to  5  because it is intended to work from subspecises to family

  # cycle to scale the traits among taxonomic levels
  for (i in 10:5) {
    temp <- ref[, -which(c(1:10) != i)]
    names(temp)[1] <- "Taxa"
    temp <- temp[temp[, 1] != "", ]
    ref.name <- ref[rownames(ref) %in% rownames(temp), 5:10]
    ref.name <- apply(ref.name, 1, function(x) (rev(x)[rev(x) != ""][1]))
    ref.na <- c(ref.na, ref.name)
    trait.interm <- rbind(trait.interm, temp)
  }

  trait_db <- trait.interm

  DF <- x[["Tree"]]

  # allow the user to work at a desired taxonomic level

  if (!tax_lev %in% c(
    "Family", "Subfamily", "Tribus", "Genus",
    "Species", "Subspecies", "Taxa"
  )) {
    stop("Maximum taxonomic level is family.")
  }

  if (!identical(tax_lev, "Taxa")) {
    # set values to "" for the taxonomic levels lower than specified
    DF <- DF[, 1:11]
    DF[colnames(DF)[(which(colnames(DF) %in% tax_lev) + 1):11]] <- ""

    # remove duplicated rows
    DF <- DF[!duplicated(DF), ]

    # change Taxa to the taxa of the required taxonomic level
    DF[, 11] <- DF[, tax_lev]
  }

  # DFtaxa stores the taxa present in the user database

  DFtaxa <- as.character(DF[11])
  result.list <- apply(DF, 1, function(x) traitS(x = x, y = DF, z = trait_db, w = ref.na))
  result.data.frame <- do.call(rbind, result.list)
  result.data.frame <- result.data.frame[result.data.frame$Taxa_db != "", ]
  result.data.frame.single <- result.data.frame

  unique.taxa <- unique(result.data.frame.single$Taxa_db)

  # deleting row for the following reason. It happens that traits are present for a species
  # and also for the genus of this species. If the user sample
  # contains a species other than that reported in the trait list we want only to keep the
  # trait at genus level

  for (i in 1:length(unique.taxa)) {
    res <- result.data.frame.single[result.data.frame.single$Taxa_db %in% unique.taxa[i], 1:3]
    res.sum <- (res[, 1] == res[, 2]) + (res[, 1] == res[, 3]) + (res[, 2] == res[, 3])
    if (any(res.sum == 0)) {
      if (sum(res.sum) != 0) {
        to.del <- which(result.data.frame.single$Taxa_db %in% unique.taxa[i])[res.sum == 0]
        result.data.frame.single <- result.data.frame.single[-to.del, ]
      }
    }
  }

  ref_long <- data.frame(Taxonomic_level = character(), Taxa = character(), stringsAsFactors = FALSE)

  for (i in 10:1) {
    temp <- as.character(dfref[, i])
    temp <- temp[temp != ""]
    temp.rep <- rep(names(ref[, i, drop = FALSE]), length(temp))
    temp.df <- data.frame(Taxonomic_level = temp.rep, Taxa = temp)
    ref_long <- rbind(ref_long, temp.df)
  }

  ref_long <- ref_long[!duplicated(ref_long), ]

  taxa_db.taxlev <- ref_long[match(ref_long$Taxa, result.data.frame.single$Taxa_db), ]

  result.data.frame.single$Taxa_db <- as.character(result.data.frame.single$Taxa_db)
  result.data.frame.single$Traits_real <- as.character(result.data.frame.single$Traits_real)
  ref_long$Taxa <- as.character(ref_long$Taxa)

  taxa_db.taxlev <- inner_join(result.data.frame.single[, 1, drop = FALSE], ref_long, by = c("Taxa_db" = "Taxa"))
  traits.taxlev <- inner_join(result.data.frame.single[, 2, drop = FALSE], ref_long, by = c("Traits_real" = "Taxa"))
  names(taxa_db.taxlev) <- c("Taxa_taxlev", "Taxa_db")
  names(traits.taxlev) <- c("Traits_taxlev", "Traits_real")


  tax_lev.info <- data.frame(taxa_db.taxlev, traits.taxlev, stringsAsFactors = FALSE)

  dist.taxlev <- data.frame(tax_lev = names(dfref)[-11], distance = c(10:1))
  dist.taxlev$tax_lev <- as.character(dist.taxlev$tax_lev)
  tax_lev.info$Taxa_db <- as.character(tax_lev.info$Taxa_db)
  tax_lev.info$Traits_real <- as.character(tax_lev.info$Traits_real)

  a <- inner_join(tax_lev.info[, c(2, 4)], dist.taxlev, by = c("Taxa_db" = "tax_lev"))[, 3]
  b <- inner_join(tax_lev.info[, c(2, 4)], dist.taxlev, by = c("Traits_real" = "tax_lev"))[, 3]

  tax_lev.info$Taxonomic_distance <- a - b

  names(tax_lev.info)[1] <- "Taxa"
  final.traits <- data.frame(tax_lev.info, result.data.frame.single[, -c(1:3)])
  rownames(final.traits) <- NULL
  if (is.null(filter_by_distance)) {
    final.traits
  } else {
    if (is.character(filter_by_distance)) {
      if (filter_by_distance == "pos") {
        final.traits[final.traits$Taxonomic_distance >= 0, ]
      } else {
        if (filter_by_distance == "neg") {
          final.traits[final.traits$Taxonomic_distance <= 0, ]
        } else {
          stop("pos, neg or an integer are needed when filter_by_distance is not NULL")
        }
      }
    } else {
      final.traits[abs(final.traits$Taxonomic_distance) <= abs(filter_by_distance), ]
    }
  }
}
