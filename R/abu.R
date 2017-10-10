abu <- function(x){

  # check if the object d is of class "biomonitoR"

  if (class(x) != "biomonitoR") {
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    return("Object x is not an object of class biomonitoR")
  }

  nord <- apply(x[["Order"]][,-1, , drop = F], 2, FUN=sum)
  return(nord)
}
