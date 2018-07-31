checkNames <- function(x){

  if(group == "mi"){
    dictio <- system.file("dict", "mi_dictionary.txt", package="biomonitoR")
  }
  if(group == "mf"){
    dictio <- system.file("dict", "mf_dictionary.txt", package="biomonitoR")
  }
  taxaCar <- as.character(x$Taxa)

  # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
  # hunspell_check and hunspell_suggest are from the package hunspell
  nameCheck <- hunspell::hunspell_check(taxaCar, dict = dictio)
  nameSuggest <- hunspell::hunspell_suggest(taxaCar, dict = dictio)
  sum(nameCheck)==length(nameCheck)
}
