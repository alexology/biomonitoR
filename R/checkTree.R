checkTree <- function(x){
  x.names <- names( x )
  
  # reverse x.names order while removing Taxa
  x.check <- rev(x.names)[-1]
  
  # check if there are duplicates in Taxa column of the tree provided by the user
  
  # initialize the error message
  
  mes <- c()
  
  if( nrow( x ) != length( unique( x[ , "Taxa"] ) ) ){
    x.taxa <- x[ , "Taxa"]
    taxa.dup <- x.taxa[duplicated(x.taxa)]
    taxa.dup <- paste(taxa.dup, collapse = ",")
    mes <- paste("Duplicate names are present in Taxa column:", taxa.dup)
  }
  
  
  row.count <- c()
  
  for(i in 1:length(x.check)){
    t.t <- x.check[ i ]
    if(is.null( row.count ) == T){
      temp.df <- x
    } else {
      temp.df <- x[ -row.count, ]
    }
    
    t.names <- temp.df[,  t.t]
    t.nrow <- which(temp.df[, t.t] != "")
    t.names <- t.names[t.nrow]
    
    # check if duplicate names are present in the columns of the tree provided by the user
    
    if((length(unique(t.names)) - length(t.names)) > 0){
      mes <- paste("Duplicate names are present in ", t.t, " column")
    }
    
    t.check <- sum(x[ t.nrow, "Taxa"] == length(t.names))
    if(t.check != 0){
      t.diff <- setdiff(t.names, x[ t.nrow, "Taxa"])
      mes <- paste("Columns Taxa and", t.t, "differ by: ", t.diff )
    }
    
    row.count <- c(row.count, t.nrow)
    
  }
  return(mes)
}