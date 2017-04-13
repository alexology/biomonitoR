checkNames <- function(x){
  dictio <- system.file("dict", "ephemeroptera_ref.txt", package="biomonitoR")
  taxaCar <- as.character(x$Taxa)

  # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
  # hunspell_check and hunspell_suggest are from the package hunspell
  nameCheck <- hunspell::hunspell_check(taxaCar, dict = dictio)
  nameSuggest <- hunspell::hunspell_suggest(taxaCar, dict = dictio)
  sum(nameCheck)==length(nameCheck)
}
