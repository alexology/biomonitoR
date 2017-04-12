checkNames <- function(x){
  dictio <- paste(find.package("biomonitoR"),"/ephemeroptera_ref.txt", sep="")
  taxaCar <- as.character(x$Taxa)

  # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
  # hunspell_check and hunspell_suggest are from the package hunspell
  nameCheck <- hunspell_check(taxaCar, dict = dictio)
  nameSuggest <- hunspell_suggest(taxaCar, dict = dictio)
  sum(nameCheck)==length(nameCheck)
}
