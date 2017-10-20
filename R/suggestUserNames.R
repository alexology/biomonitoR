suggestUserNames <- function(x){
  dictio <- system.file("dict", "macro_dictionary.txt", package="biomonitoR")
  taxaCar <- as.character(x$Taxa)
  taxaCar <- as.character(x$Taxa)
  
  # nameCheck and nameSuggest check for the wrong names and suggest for correct names.
  # hunspell_check and hunspell_suggest are from the package hunspell
  nameCheck <- hunspell::hunspell_check(taxaCar, dict = dictio)
  nameSuggest <- hunspell::hunspell_suggest(taxaCar, dict = dictio)
  
  
  n <- which(nameCheck==F) #number of wrong names
  if(length(n) == 0){
    return(n)
  }
  else{
  temp <- rep(NA, length(n)) # vector to store user choices
  wrongName <- as.vector(x[n,"Taxa"]) # vector of wrong taxa names
  correctName <- nameSuggest[n] # list of suggestions
  for(i in 1:length(n)){
    choice <- c(correctName[[i]],"User taxon input") # choices provided to the user
    # provide suggestions to the user and if the user can't find the correct name he can enter
    # the right name by himself
    temp[i] <- select.list(choice, title=wrongName[i])
    if(temp[i]=="User taxon input"){
      temp[i] <- readline(prompt="Enter taxon name: ")
    }
    else{
      temp[i] <- temp[i]
    }
  }
  # Remove hashtag for a standalone function
  return(list("wrongNames" = wrongName, "correctNames" = temp))
  }
}
