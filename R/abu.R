abu <- function(x) {

  # calculate abundance from the Taxa level. At this level
  # there should not be loss of information. The -1 is to remove the
  # taxa column. Two option because of the presence-absence issue

  res <- apply(x[["Taxa"]][, -1, drop = FALSE], 2, FUN = sum)

  res
}
