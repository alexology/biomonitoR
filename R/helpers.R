om <- function(x) {
  ifelse(is.null(x), NA, x)
}

# Only characters
only_char <- function(x) {
  ifelse(grepl("[^A-Za-z]", gsub(" ", "", x)), NA, x)
}


# Count words
countWords <- function(x){
  length(strsplit(x,' ')[[1]])
}

# Species name
specieWorms <- function(x){
  if(countWords(only_char(om(x))) == 1){
    return(NA)
  }

  if(countWords(only_char(om(x))) == 2){
    return(x)
  }

  if(countWords(only_char(om(x))) == 3){
    sp <- strsplit(only_char(om(x)), " ")[[1]]
    return(paste(sp[1], sp[2]))
  }
}




# Species name
specieName <- function(x){
  if(countWords(only_char(om(x))) == 1){
    return(NA)
  }

  if(countWords(only_char(om(x))) == 2){
    return(x)
  }

  if(countWords(only_char(om(x))) == 3){
    sp <- strsplit(only_char(om(x)), " ")[[1]]
    return(paste(sp[1], sp[2]))
  }
}

# First letter upper case
upperFirst <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
