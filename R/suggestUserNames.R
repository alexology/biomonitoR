#' @importFrom utils select.list
#' @importFrom hunspell hunspell_check hunspell_suggest

suggestUserNames <- function(x, group = "mi") {
  if (group == "mi") {
    dic.path <- system.file("dict", "mi_dictionary.txt", package = "biomonitoR")
    dictio <- dictionary(dic.path, cache = FALSE)
  }
  if (group == "mf") {
    dic.path <- system.file("dict", "mf_dictionary.txt", package = "biomonitoR")
    dictio <- dictionary(dic.path, cache = FALSE)
  }
  if (group == "fi") {
    dic.path <- system.file("dict", "fi_dictionary.txt", package = "biomonitoR")
    dictio <- dictionary(dic.path, cache = FALSE)
  }
  if (group == "di") {
    dic.path <- system.file("dict", "di_dictionary.txt", package = "biomonitoR")
    dictio <- dictionary(dic.path, cache = FALSE)
  }
  taxaCar <- as.character(x$Taxa)
  taxaCar <- sapply(taxaCar, capWords, USE.NAMES = FALSE)

  # replace space with underscore to be compatible with hunspell
  taxaCar <- gsub(" ", "_", taxaCar)

  # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
  # hunspell_check and hunspell_suggest are from the package hunspell
  nameCheck <- hunspell::hunspell_check(taxaCar, dict = dictio)
  nameSuggest <- hunspell::hunspell_suggest(taxaCar, dict = dictio)

  # This part of the function change the wrong names to correct
  # the user is provided with an interactive selection interface with select.list

  n <- which(nameCheck == FALSE) # number of wrong names
  if (length(n) == 0) {
    return(n)
  }
  else {
    temp <- rep(NA, length(n)) # vector to store user choices
    wrongName <- as.vector(x[n, "Taxa"]) # vector of wrong taxa names
    correctName <- nameSuggest[n] # list of suggestions
    for (i in 1:length(n)) {
      choice <- c(correctName[[i]], wrongName[i], "User taxon input", "Delete taxon", "exit") # choices provided to the user
      # provide suggestions to the user and if the user can't find the correct name he can enter
      # the right name by himself
      temp[i] <- select.list(choice, title = wrongName[i])
      if (temp[i] == "User taxon input") {
        temp[i] <- readline(prompt = "Enter taxon name: ")
      }
      else {
        if (temp[i] == "Delete taxon") {
          temp[i] <- "REMOVe"
        } else {
          if (temp[i] == "exit") {
            stop()
          } else {
            temp[i] <- temp[i]
          }
        }
      }
    }

    return(list("wrongNames" = wrongName, "correctNames" = temp))
  }
}
