#' refFromTree
#'
#' @description
#' \Sexpr[results=rd, stage=render]{ lifecycle::badge("deprecated") }
#'
#' This function transforms a taxonomic tree to a reference database suitable for biomonitoR.
#' @param x taxonomic tree. See \code{\link{Tree}} for an example.
#' @param group merge the user database with one of the biomonitoR reference databases, default to `none`.
#'  If duplicated Taxa names are present this function keeps the name provided by the user.
#'  Check the reliability of results when using group = `mi` for macroinvertebrates or group = `mf` for macrophytes.
#' @keywords refFromTree
#' @export


refFromTree <- function(x, group = "none") {
  .Deprecated("ref_from_tree", package = "biomonitoR")

  n <- ncol(x)

  if (n == 1) {
    stop("data.frame with 1 column are not allowed!")
  }

  # dummy variables to avoid RCMD notes

  Phylum <- Class <- Subclass <- Order <- Family <- Subfamily <- Tribus <- Genus <- Species <- Subspecies <- NULL

  taxa.col <- c(
    "Phylum", "Class", "Subclass", "Order", "Family", "Subfamily",
    "Tribus", "Genus", "Species", "Subspecies", "Taxa"
  )

  empty.df <- as.data.frame(matrix(ncol = length(taxa.col), nrow = 0))
  colnames(empty.df) <- taxa.col
  DF <- empty.df



  # transform all the columns to character and remove leading and trailing whitespaces
  x[] <- lapply(x, as.character)
  x[] <- lapply(x, trimws)

  # remove duplicated rows
  x <- x[!duplicated(x), ]

  # colnames of the reference database provided by the user
  cref.name <- colnames(x)

  # transform to capital letter
  cref.name <- sapply(cref.name, capWords, USE.NAMES = FALSE)

  # test if cref.name are different from colnames accepted by biomonitoR
  if (sum(!cref.name %in% taxa.col) != 0) {
    stop("Provide valid column names")
  }

  # reorder user data.frame according to biomonitoR taxa tree
  value.match <- match(taxa.col, names(x))
  x <- x[, value.match[!is.na(value.match)]]

  cref.name <- colnames(x)

  for (i in 1:n) {
    temp.name <- colnames(x[, i, drop = FALSE])
    temp.pos <- which(cref.name == temp.name)
    temp <- x [, 1:temp.pos, drop = FALSE]

    # remove rows with empty cells
    temp <- temp[which(temp[, temp.name] != ""), , drop = FALSE]
    temp.un <- unique(temp)
    empty.df <- merge(empty.df, temp.un, all.y = TRUE, sort = FALSE)
    empty.df$Taxa <- temp.un[, temp.name]
    DF <- rbind.data.frame(DF, empty.df)
  }

  # remove NA and transform characters to factors with as.data.frame
  DF[is.na(DF)] <- ""
  DF <- as.data.frame(unclass(DF))

  # remove leading and final spaces
  DF <- sapply(DF, trim, USE.NAMES = FALSE)

  # remove NA originating from empty columns when doing as.data.frame
  if (length(DF[is.na(DF)]) > 0) {
    DF[is.na(DF)] <- ""
  }

  # check for duplicates or errors

  DF <- as.data.frame(DF)
  DF <- DF[!duplicated(DF), ]

  s.mes <- checkTree(DF)
  if (is.null(s.mes) == FALSE) {
    stop(s.mes)
  }

  if (group == "mi") {
    DF <- rbind(DF, mi_ref)
    DF <- DF[!duplicated(DF$Taxa), ]
  }

  if (group == "mf") {
    DF <- rbind(DF, mf_ref)
    DF <- DF[!duplicated(DF$Taxa), ]
  }

  DF <- DF[, match(taxa.col, colnames(DF))]
  DF
}
