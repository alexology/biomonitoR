#' @describeIn allindices Simpson's Index of Diversity

simpson <- function(x, tax_lev = "Taxa") {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  # get the data.frame at the desired taxonomic level
  DF <- x[[tax_lev]]

  if (inherits(x, "bin")) {
    DF <- to_bin(DF)
  }

  # remove unassigned row from the species count if present
  if ("unassigned" %in% DF[, 1]) {
    z <- which(DF[, 1] == "unassigned")
    DF <- DF[-z, ]
  }

  # apply the function for calculating the Simpson's Index of Diversity that is stored in the Pi.R file
  res <- apply(DF[, -1, drop = FALSE], 2, FUN = function(x) {
    Pi(x, index = "Simpson")
  })
  res
}
