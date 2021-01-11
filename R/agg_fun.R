# helper function for combTaxa
# it takes as input a list of which each element is a subset of the data.frame with the x-th taxa combination
# the first column of the input is thus composed by a vector of taxa (e.g. Simuliidae and Tanypodinae if n is equal to 2) while the other columns
# store the abundance of the taxa for all the samples in the data.frame
# this function sum the abundances of the x-th combination of and change the label merging the taxa names (e.g. Simuliidae_Tanypodinae)

agg_fun <- function(x, rel = FALSE) {
  lab <- paste(x[, 1], collapse = "_")
  temp <- apply(x[, -1], 2, FUN = sum)
  DF <- data.frame(Taxa = lab, t(temp))
  return(DF)
}
