#' @importFrom stats as.dist

ddis <- function(x, st.names) {
  com <- x[, -which(names(x) %in% st.names)]
  n <- nrow(com)
  diff.ar.ij <- function(i, j, data) {
    sum(data[, i] != data[, j])
  }
  corp <- Vectorize(diff.ar.ij, vectorize.args = list("i", "j"))
  temp <- outer(1:n, 1:n, corp, data = t(com))
  temp <- temp / max(temp) * 100
  rownames(temp) <- rownames(com)
  return(temp)
}

delta <- function(x, dis) {
  a <- as.dist(outer(x, x))
  b <- a * as.dist(dis)
  d <- b / sum(x) / (sum(x) - 1) * 2
  return(sum(d))
}

delta.st <- function(x, dis) {
  a <- as.dist(outer(x, x))
  b <- a * as.dist(dis)
  d <- b / sum(a)
  return(sum(d))
}
