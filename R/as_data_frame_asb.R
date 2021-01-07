#' @export
as.data.frame.asb <- function(x, row.names, optional, ..., object = 1, stringsAsFactors = FALSE) {
  x.name <- deparse(substitute(x))

  if (sum(object %in% seq_along(x)) == 0) stop(paste(x.name, " is of length ", length(x), sep = ""))
  DF <- unclass(x[[object]])
  i <- sapply(DF, is.factor)
  DF[i] <- lapply(DF[i], as.character)
  as.data.frame(DF, stringsAsFactors = stringsAsFactors, ...)
}


#' @export
print.asb <- function(x, ...) {
  DF <- unclass(x)
  if (length(DF) == 1) print(as.data.frame(x))
  if (length(DF) > 1) print(DF)
}


#' subset.asb
#'
#' @description `summary` method for class `asb`.
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("experimental") }
#'
#' @param x results of function `as_biomonitor()`.
#' @param ... further arguments to be passed to or from other methods.
#' @param taxa a taxon or a vector of taxa to retain.
#' @param exclude a taxon or a vector of taxa to exclude. For example, this option is useful to exclude alien taxa.
#'
#' @export
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#'
#' # select EPT Taxa (Ephemeroptera, Plecoptera and Trichoptera)
#'
#' subset(data_bio, taxa = c("Ephemeroptera", "Plecoptera", "Trichoptera"))
#'
#' # select Trichoptera excluding the trichopteran family Hydropsychidae
#'
#' tricho <- subset(data_bio, taxa = "Trichoptera", exclude = "Hydropsychidae")
#' tricho_agg <- aggregate_taxa(tricho)
#'
#' # exclude Chironomidae
#' subset(data_bio, exclude = "Chironomidae")


subset.asb <- function(x, ..., taxa = NULL, exclude = NULL) {
  DF <- as.data.frame(x)

  if (is.null(taxa) & is.null(exclude)) {
    stop("taxa and exclude cannot be both null")
  }

  # # check if the taxa argument is empty or contains NULL strings
  # if ( is.null( taxa ) == TRUE || ( any( taxa == "" ) & length( taxa ) == 1 ) ) {
  #    stop( "Please provide at least taxon name" )
  # }

  # transform the taxa argument to character

  if (!is.null(taxa)) taxa <- trimws(as.character(taxa))
  if (!is.null(exclude)) exclude <- trimws(as.character(exclude))

  # get the list of all the taxa present in the x object including the information stored
  # in the taxonomic tree

  df.vec <- unique(as.character(unlist(DF[, 1:11])))
  df.vec <- df.vec[!df.vec %in% ""]

  # get the name of the x object
  x.name <- deparse(substitute(x))

  # find any taxa not present in the dataset provided by the user. If no taxa belongs to the list of all the taxa present
  # in the x object the function will be stopped, otherwise it will provide the list of missing taxa
  if (!is.null(taxa) & is.null(exclude)) {
    taxa.sub <- taxa[!taxa %in% df.vec]
  }

  if (is.null(taxa) & !is.null(exclude)) {
    taxa.sub <- exclude[!exclude %in% df.vec]
  }

  if (!is.null(taxa) & !is.null(exclude)) {
    comb.taxa <- c(taxa, exclude)
    taxa.sub <- comb.taxa[!comb.taxa %in% df.vec]
  }

  if (length(taxa.sub) > 0) {
    print(paste("The following taxon were not find in the ", x.name, " database and has been excluded: ", taxa.sub, sep = ""))

    if (!is.null(taxa) & is.null(exclude)) {
      taxa <- taxa[taxa %in% df.vec]
    }

    if (is.null(taxa) & !is.null(exclude)) {
      exclude <- exclude[exclude %in% df.vec]
    }

    if (!is.null(taxa) & !is.null(exclude)) {
      taxa <- taxa[taxa %in% df.vec]
      exclude <- exclude[exclude %in% df.vec]
    }

    if (length(taxa) == 0 & length(exlcude) == 0) {
      stop("None of the taxa provided were found in the ", x.name, " database")
    }
  }

  # list the row and column numbers at which the desired taxa are stored
  if (!is.null(taxa)) {
    taxind <- data.frame(row = numeric(), col = numeric())
    for (i in seq_along(taxa)) {
      temp <- which(DF == taxa[i], arr.ind = T)
      taxind <- rbind(temp, taxind)
    }

    to.keep <- unique(taxind$row)
  }


  if (!is.null(exclude)) {
    # list the row and column numbers at which the desired taxa are stored
    taxind.e <- data.frame(row = numeric(), col = numeric())
    for (i in seq_along(exclude)) {
      temp <- which(DF == exclude[i], arr.ind = T)
      taxind.e <- rbind(temp, taxind.e)
    }
    to.rem <- unique(taxind.e$row)
  }

  if (!is.null(taxa) & !is.null(exclude)) {
    to.keep <- to.keep
  }

  if (is.null(taxa) & !is.null(exclude)) {
    to.keep <- setdiff(1:nrow(DF), to.rem)
  }

  if (!is.null(taxa) & !is.null(exclude)) {
    to.keep <- setdiff(to.keep, to.rem)
  }


  df_sub <- list(DF[to.keep, ])
  names(df_sub) <- ""
  class(df_sub) <- class(x)
  df_sub
}
