om <- function(x) {
  ifelse(is.null(x), NA, x)
}

# Only characters
only_char <- function(x) {
  ifelse(grepl("[^A-Za-z]", gsub(" ", "", x)), NA, x)
}
