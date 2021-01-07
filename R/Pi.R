#' @importFrom stats optimize

Pi <- function(x, index = "Shannon", base = exp(1)) {
  ntaxa <- length(x[x > 0]) # calculates the taxa rixhness
  abu <- sum(x) # calculates total abundance of the vector
  relAbu <- x / abu # calculates the relative abundance
  if (index == "Shannon") {
    temp <- relAbu * log(relAbu, base)
    indValue <- -sum(temp[is.finite(temp)])
  }
  if (index == "Simpson") {
    indValue <- 1 - sum(relAbu^2)
  }
  if (index == "Simpsoneven") {
    indValue <- (1 / sum(relAbu^2)) / ntaxa
  }
  if (index == "Invsimpson") {
    indValue <- 1 / sum(relAbu^2)
  }
  if (index == "Margalef") {
    indValue <- (ntaxa - 1) / log(abu)
  }
  if (index == "Menhinick") {
    indValue <- ntaxa / sqrt(abu)
  }
  if (index == "Berpar") {
    indValue <- max(x) / abu
  }
  if (index == "Invberpar") {
    indValue <- 1 / (max(x) / abu)
  }
  if (index == "Brillouin") {
    indValue <- (lfactorial(abu) - sum(lfactorial(x))) / abu
  }
  if (index == "Mcintosh") {
    indValue <- (abu - sqrt(sum(x^2))) / (abu - sqrt(abu))
  }
  if (index == "Fisher") {
    fn <- function(r, n, x) (abs((1 - x) * n / x * (-log(1 - x)) - r))
    temp <- optimize(f = fn, interval = c(0, 1), n = abu, r = ntaxa, tol = 0.000000001)[[1]]
    indValue <- abu * (1 - temp) / temp
  }
  indValue
}
