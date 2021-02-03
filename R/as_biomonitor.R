#' as_biomonitor
#'
#' @description
#' This function merge the user dataframe with a reference database and suggest corrections for mispelled names.
#'
#' @param x a data.frame with a column called "Taxa" where store taxa names and samples on the other columns (see the example macro_ex).
#' @param group biotic group of interest. Possible values are `mi` for macroinvertebrates, `mf` for macrophytes and `fi` for fish. The choice will set the right reference database for the specified group.
#'  This option will not be considered if a custom reference database is provided. Default to `mi`.
#' @param dfref a custom reference database that replaces the reference database.
#' @param to_change a `data.frame` specifying the taxa name that needs to be changed.
#' This `data.frame` needs a column called *Taxon* containing the taxon to aggregate and a column called *Correct_Taxon* with the aggregation specifications.
#' By default, when group is set to `mi` Hydracarina, Hydracnidia and Acariformes are changed to Trombidiformes.
#' @param FUN the function to be applied for aggregating rows with duplicated taxa names.
#' It should be `sum` for abundances, while it should be `bin` for presence-absence data. Default to `sum`.
#' @param correct_names if `TRUE` alternative names will be suggested for taxa not found in the reference database. Default to
#' `FALSE`, with which the unrecognized taxa will be removed.
#' @param traceB track changes in taxa names.
#'
#' @keywords as_biomonitor
#' @details The function `as_biomonitor()` checks the taxonomy of the data.frame provided by the user and suggests correction for mispelled names.
#' If one or more taxa names of `x` are not present in the reference database or the spell checker is not able to find any suggestion the user is asked to exit.
#' This behaviour is to assure consistency with other functions implemented in biomonitoR.
#' Default references databases are provided for macroinvertebrates and macrophytes.
#' Both databases heavily rely on the information provided by the [freshwaterecology.info](https://www.freshwaterecology.info/) website.
#' If `dfref` is not NULL a custom dictionary will be saved in the working directory to let the `as_biomonitor ()` function work correctly.
#' If you are unable to build a reference database by your own please check the function \code{\link{ref_from_tree}} for a possible solution.
#' `as_biomonitor()`  returns an object of class `asb` togheter with one of the classes `abundance` or `bin`.
#' The function \code{\link{quick_rename}} works as the `as_biomonitor()` but returns
#' a data.frame without the biomonitoR format.
#' `as_biomonitor()` aggregates all the rows with the same name with the option `FUN` and converts all the `NA` to 0.
#' If only 1 and 0 are present `x` will be imported as presence-absence.
#' When `group = mi` Hydracarina, Hydracnidia or Acariformes are changed to Trombidiformes given the uncertain taxonomic status of this group.
#'
#' @importFrom stats aggregate
#' @importFrom utils select.list stack
#' @importFrom hunspell dictionary hunspell_check hunspell_suggest
#' @export
#' @seealso [quick_rename] [ref_from_tree]
#' @references Schmidt-Kloiber, A., & Hering, D. (2015). www.freshwaterecology.info -
#' An online tool that unifies, standardises and codifies more than
#' 20,000 European freshwater organisms and their ecological preferences.
#' Ecological indicators, 53, 271-282.
#' @examples
#' data(mi_prin)
#' data_bio <- as_biomonitor(mi_prin, group = "mi")
as_biomonitor <- function(x, group = "mi", dfref = NULL, to_change = "default", FUN = sum, correct_names = FALSE, traceB = FALSE) {

  # check if user database contains a column called Taxa
  if (!"Taxa" %in% names(x)) {
    stop("A column called Taxa is needed")
  }

  asb.call <- as.character(as.list(match.call())[["FUN"]])
  if (length(asb.call) == 0) {
    asb.call <- "sum"
  }

  # check if columns other than Taxa are numeric
  # position of column Taxa
  col.taxa <- which(names(x) == "Taxa")
  col.class <- sapply(x[, -col.taxa], is.numeric)
  if (any(col.class == FALSE)) {
    stop("Non-numeric columns other than Taxa are not allowed")
  }

  # check if NAs are present. If present they are changed to 0.
  check.na <- any(is.na(x))
  x[is.na(x)] <- 0

  # initialize check.pa
  check.pa <- FALSE

  if (!is.null(to_change) & !is.data.frame(to_change) & !identical(to_change, "default")) (stop("to_change needs to be NULL, default or data.frame as specified in the help"))

  if (is.data.frame(to_change)) {
    if (length(unique(to_change$Taxon)) != nrow(to_change)) (stop("the same name cannot be present twice in the Taxon column of the data.frame to_change"))
    if (nrow(to_change) == 0) (stop("to_change must have at least one entry"))
  }

  if (!any(x[, !colnames(x) %in% "Taxa"] > 1) & all(x[, !colnames(x) %in% "Taxa"] %% 1 == 0) & !identical(asb.call, "bin")) {
    (warning("Presence-absence data detected but FUN is not set to bin. Is it this what you want?"))
    check.pa <- TRUE
  }


  if (any(x[, !colnames(x) %in% "Taxa"] > 1) & all(x[, !colnames(x) %in% "Taxa"] %% 1 == 0) & !identical(asb.call, "sum")) (warning("Abundance data detected but FUN is not set to sum. Is it this what you want?"))
  if (any(x[, !colnames(x) %in% "Taxa"] %% 1 != 0)) warning("Decimal numbers detected. Please check carefully which FUN to use.")






  # set the reference database for the specified group
  if (group == "mi") {
    # macroinvertebrates
    ref <- mi_ref
  }

  if (group == "mf") {
    # macrophytes
    ref <- mf_ref
  }

  if (group == "fi") {
    # fish
    ref <- fi_ref
  }

  # allow the users to use their own reference database
  if (!is.null(dfref)) {
    dfref[is.na(dfref)] <- ""
    dfref <- as.data.frame(unclass(dfref))
    ref <- dfref
    newDictio(ref)
    group <- "custom"
  }

  # select or create dictinoaries

  if (!identical(group, "custom")) {
    if (identical(group, "mi")) {
      dic.path <- system.file("dict", "mi_dictionary.txt", package = "biomonitoR")
      # very important to set cache equal to FALSE, otherwise suggestNames will provide inconsistent results.
      dictio <- dictionary(dic.path, cache = F)
    }
    if (identical(group, "mf")) {
      dic.path <- system.file("dict", "mf_dictionary.txt", package = "biomonitoR")
      dictio <- dictionary(dic.path, cache = F)
    }
    if (identical(group, "fi")) {
      dic.path <- system.file("dict", "fi_dictionary.txt", package = "biomonitoR")
      dictio <- dictionary(dic.path, cache = F)
    }
  }

  if (identical(group, "custom")) {
    dic.path <- c(paste(getwd(), "/custom_dictio.dic", sep = ""))
    dictio <- dictionary(dic.path, cache = F)
  }


  # change Taxa from factor to character
  x$Taxa <- as.character(x$Taxa)

  # change the name of taxa to lowercase and capital letter
  x$Taxa <- trimws(sapply(x$Taxa, capWords, USE.NAMES = FALSE))


  # change the Hydracarina, Hydracnidia or Acariformes changed to Trombidiformes
  if (identical(to_change, "default")) {
    to_change_mi[, "Taxon"] <- trimws(sapply(to_change_mi[, "Taxon"], capWords, USE.NAMES = FALSE))
    to_change_mi[, "Correct_Taxon"] <- trimws(sapply(to_change_mi[, "Correct_Taxon"], capWords, USE.NAMES = FALSE))

    if (any(to_change_mi[, "Taxon"] %in% x$Taxa)) {

      # store results for traceB
      to_store <- to_change_mi[to_change_mi$Taxon %in% x$Taxa, ]

      change_uni <- x[x$Taxa %in% to_change_mi$Taxon, "Taxa", drop = TRUE]

      for (i in 1:length(change_uni)) {
        x[x$Taxa %in% change_uni[i], "Taxa"] <- to_change_mi[to_change_mi$Taxon %in% change_uni[i], "Correct_Taxon"]
      }
    }
  }

  # change according to user needs
  if (is.data.frame(to_change)) {
    to_change[, "Taxon"] <- trimws(sapply(to_change[, "Taxon"], capWords, USE.NAMES = FALSE))
    to_change[, "Correct_Taxon"] <- trimws(sapply(to_change[, "Correct_Taxon"], capWords, USE.NAMES = FALSE))

    if (any(to_change[, "Taxon"] %in% x$Taxa)) {

      # store results for traceB
      to_store <- to_change[to_change$Taxon %in% x$Taxa, ]


      change_uni <- x[x$Taxa %in% to_change$Taxon, "Taxa", drop = TRUE]

      for (i in 1:length(change_uni)) {
        x[x$Taxa %in% change_uni[i], "Taxa"] <- to_change[to_change$Taxon %in% change_uni[i], "Correct_Taxon"]
      }

    }
  }


  # search for mispelled names
  ref_taxa <- unique(ref$Taxa)
  wrong_taxa <- unique(x$Taxa)

  # get wrong names
  wrong_taxa <- wrong_taxa[!wrong_taxa %in% ref_taxa]

  if (length(wrong_taxa) > 0) {
    # replace space with underscore to be compatible with hunspell
    wrong_taxa <- gsub(" ", "_", wrong_taxa)

    # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
    # hunspell_check and hunspell_suggest are from the package hunspell

    name_suggest <- hunspell_suggest(wrong_taxa, dict = dictio)

    names(name_suggest) <- wrong_taxa

    names_suggest <- stack(name_suggest)
    names_suggest$ind <- as.character(names_suggest$ind)
    colnames(names_suggest) <- c("suggested", "excluded")
    names_suggest <- names_suggest[, 2:1]
  }


  ### INTERACTIVE EQUAL TO FALSE

  # if correct names is equal to FALSE all the wrong names will be discarded

  if (!correct_names) {
    x <- x[!x$Taxa %in% wrong_taxa, , drop = FALSE]
  } else {
    if (length(wrong_taxa) == 0) {
      x <- x
    } else {
      temp <- rep(NA, length(wrong_taxa)) # vector to store user choices

      for (i in 1:length(wrong_taxa)) {
        choice <- names_suggest[names_suggest$excluded %in% wrong_taxa[i], 2, drop = TRUE] # choices provided to the user
        # provide alternative names to the user and if the user can't find the correct name he must exit
        temp[i] <- select.list(choice, title = wrong_taxa[i])
        if (temp[i] == "exit") stop()
      }

      taxa_corrected <- data.frame("wrong_names" = wrong_taxa, "correct_names" = temp, stringsAsFactors = FALSE)

      # change wrong names
      for (i in 1:nrow(taxa_corrected)) {
        x[x$Taxa %in% taxa_corrected$wrong_names[i], "Taxa"] <- taxa_corrected$correct_names[i]
      }
    }
  }

  x <- aggregate(. ~ Taxa, data = x, FUN = FUN)


  if (check.pa) {
    x <- data.frame(x[, 1, drop = FALSE], (x[, -1, drop = FALSE] > 0) * 1, stringsAsFactors = FALSE)
    message("data imported as presence absence")
  }

  if (check.na) {
    message("NA detected, transformed to 0")
  }

  # merge reference database to the user data.frame
  taxa_def <- merge(ref, x, by = "Taxa", all = FALSE)

  # reorder the columns
  taxa_def <- taxa_def[, c(2:11, 1, 12:ncol(taxa_def)), drop = FALSE]



  if (!traceB) {
    taxa_def <- list(taxa_db = taxa_def)
  } else {
    if (!correct_names) {
      if (length(wrong_taxa) == 0) {
        if (!exists("to_store", inherits = FALSE)) {
          taxa_def <- list(taxa_db = taxa_def)
        } else {
          taxa_def <- list(taxa_db = taxa_def, corrected_names = to_store)
        }
      } else {
        if (!exists("to_store", inherits = FALSE)) {
          taxa_def <- list(taxa_db = taxa_def, suggested_taxa_names = names_suggest)
        } else {
          names(to_store) <- c("wrong_names", "correct_names")
          rownames(to_store) <- NULL
          taxa_def <- list(taxa_db = taxa_def, corrected_names = to_store, suggested_taxa_names = names_suggest)
        }
      }
    } else {
      if (length(wrong_taxa) == 0) {
        if (!exists("to_store", inherits = FALSE)) {
          taxa_def <- list(taxa_db = taxa_def)
        } else {
          taxa_def <- list(taxa_db = taxa_def, corrected_names = to_store)
        }
      } else {
        if (!exists("to_store", inherits = FALSE)) {
          taxa_def <- list(taxa_db = taxa_def, suggested_taxa_names = names_suggest)
        } else {
          names(to_store) <- c("wrong_names", "correct_names")
          taxa_corrected <- rbind(to_store, taxa_corrected)
          rownames(to_store) <- NULL
          taxa_def <- list(taxa_db = taxa_def, corrected_names = taxa_corrected)
        }
      }
    }
  }


  class(taxa_def) <- c("asb")

  if (identical(asb.call, "bin")) {
    class(taxa_def) <- c(class(taxa_def), "bin")
  }

  if (!is.null(dfref)) {
    class(taxa_def) <- c(class(taxa_def), "custom")
  }

  if( length(wrong_taxa) > 0 & ! traceB){
    message("Some taxa were excluded, check with traceB = TRUE")
  }

  taxa_def
}
