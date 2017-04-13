capWords <- function(x) {
  # see ?tolower
  paste(toupper(substring(x, 1,1)), substring(x, 2),
        sep="", collapse=" ")
}
