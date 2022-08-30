bray_curtis <- function(x){
  com <- x[, -unlist(lapply(x, is.numeric))]
  n <- nrow(com)
  diff.ar.ij <- function(i, j, data) {
    1 - 2*sum(pmin(data[i,], data[j, ])) / (sum(data[i, ]) + sum(data[j, ]))
  }
  corp <- Vectorize(diff.ar.ij, vectorize.args = list("i", "j"))
  temp <- outer(1:n, 1:n, corp, data = (com))
  rownames(temp) <- colnames(temp) <- rownames(com)
  as.dist(temp)
}
