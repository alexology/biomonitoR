capWords <- function(x) {
  # see ?tolower
  cap <- paste(toupper(substring(x, 1,1)), tolower(substring(x, 2)),
        sep="", collapse=" ")
}
