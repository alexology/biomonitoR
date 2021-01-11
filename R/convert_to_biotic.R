#' @describeIn convert_to_vegan convert data to biotic format

convert_to_biotic <- function(x) {
  classCheck(x)

  temp <- x$Tree

  # Excluding Taxa from the search, since it could bias further results
  tree <- temp[, 1:10]
  oligo.pos <- which(tree == "Oligochaeta", arr.ind = TRUE)


  if (nrow(oligo.pos) == 0) {
    DF <- x$Family
    DF$Family <- as.character(DF$Family)
  }

  oligo.sub <- tree[oligo.pos[, 1], 1:10]
  oligo.taxlev <- names(temp)[unique(oligo.pos[, 2])]

  if (nrow(oligo.pos) > 0 & sum(oligo.sub[, "Family"] != "") < nrow(oligo.sub)) {
    rem <- as.character(oligo.sub[, "Family"][oligo.sub[, "Family"] != ""])
    temp.fam <- x[["Family"]][!x[["Family"]]$Family %in% c(rem, "Oligochaeta"), ]
    temp.fam$Family <- as.character(temp.fam$Family)
    oligo.row <- x[[oligo.taxlev]][x[[oligo.taxlev]][, 1] == "Oligochaeta", ]
    names(oligo.row)[1] <- "Family"
    oligo.row[, 1] <- as.character(oligo.row[, 1])
    DF <- rbind(temp.fam, oligo.row)
  }

  if (nrow(oligo.pos) > 0 & sum(oligo.sub[, "Family"] != "") == nrow(oligo.sub)) {
    DF <- x$Family
  }

  if (nrow(oligo.pos) > 0 & sum(oligo.sub[, "Family"] == "") == nrow(oligo.sub)) {
    temp.fam <- x$Family
    temp.fam$Family <- as.character(temp.fam$Family)
    oligo.row <- x[[oligo.taxlev]][x[[oligo.taxlev]][, 1] == "Oligochaeta", ]
    names(oligo.row)[1] <- "Family"
    oligo.row[, 1] <- as.character(oligo.row[, 1])
    DF <- rbind(temp.fam, oligo.row)
  }

  if (inherits(x, "bin")) {
    DF <- to_bin(DF)
  }

  if ("unassigned" %in% DF[, 1]) {
    z <- which(DF[, 1] == "unassigned")
    DF <- DF[-z, ] # remove unassigned row from the species count
  }

  DF <- DF[order(DF[, 1]), ]
  DF[, 1] <- as.factor(DF[, 1])
  names(DF)[1] <- "Taxon"
  DF
}
