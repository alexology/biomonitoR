#' @name allindices
#' @title Functions for calculating diversity, eveness and dominance indices.
#'
#' @description Functions for calculating shannon, simpson, margalef, menhinick, pielou and other indices.
#'
#' @aliases  shannon simpson margalef menhinick pielou berpar brill esimpson invberpar invsimpson mcintosh fisher allindices
#'
#' @param x Result of the function `aggregate_taxa()`.
#' @param base The base of the logarithm
#' @param tax_lev Taxonomic level on which the calculation has to be made.
#'
#' @details Shannon index:
#' \deqn{H'=-\sum_{i=1}^{S} p_{i}\log(p_{i})}
#' Pielou index:
#' \deqn{J'=\frac{H'}{\log(S)}}
#' Margalef diversity:
#' \deqn{D_{mg}=\frac{(S - 1)}{\log(N)}}
#' Menhinick diversity:
#' \deqn{D_{mn}=\frac{S}{\sqrt(N)}}
#' Brillouin index:
#' \deqn{HB=\frac{\log(N) -\sum(\log(n_{i})) }{N}}
#' Simpson's index is calculated as D. biomonitoR returns 1 - D and 1/D when `index` is set to `Simpson` or `Invsimpson`.
#' \deqn{D=\sum_{i=1}^{S} p_{i}^{2}}
#' Simpson's evenness:
#' \deqn{D_{1/D}=\frac{1/D}{\sqrt(S)}}
#' Berger-Parker index:
#' \deqn{d=\frac{N_{max}}{N}}
#' The inverse of this index is also provided when `index` is set to `Invberpar`. McIntosh's diversity:
#' \deqn{D=\frac{N-U}{N-\sqrt(N)}}
#' where
#' \deqn{U=\sqrt(\sum_{i=1}^{S} n_{i}^{2})}
#' Fisher alpha is calculated as follow:
#' \deqn{\alpha=\frac{N(1-x)}{x}}
#' where x is estimated from the iterative solution of:
#' \deqn{\frac{S}{N}=-\frac{N(1-x)}{x} \log(1-x)}
#'
#'
#' p_i is the proportion of individuals found in the i-th species, n_i is the number of individuals found in the i-th species, S the species richness,
#'  N the number of individuals and N_max the number of individuals in the most abundant species. All the indices are calculated according to Magurran (2004).
#' @keywords shannon, simpson, margalef, menhinick, pielou
#' @export
#' @seealso [aggregate_taxa]
#' @references Magurran, A. E. (2004). Measuring biological diversity. Blackwell Science ltd.
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' allindices(data_agr)
#' shannon(data_agr)
#'
#' # base 2
#' shannon(data_agr, base = 2)
#'
#'
#' @export simpson
#' @export esimpson
#' @export invsimpson
#' @export margalef
#' @export menhinick
#' @export pielou
#' @export berpar
#' @export invberpar
#' @export brill
#' @export mcintosh
#' @export shannon
#' @export fisher

allindices <- function(x, tax_lev = "Taxa", base = exp(1)) {

  # check if the object x is of class "biomonitoR"
  classCheck(x)

  st.names <- names(x[["Taxa"]][, -1, drop = FALSE])

  f1 <- function(x) (shannon(x, tax_lev = tax_lev, base = base))
  f2 <- function(x) (berpar(x, tax_lev = tax_lev))
  f3 <- function(x) (brill(x, tax_lev = tax_lev))
  f4 <- function(x) (invberpar(x, tax_lev = tax_lev))
  f5 <- function(x) (invsimpson(x, tax_lev = tax_lev))
  f6 <- function(x) (margalef(x, tax_lev = tax_lev))
  f7 <- function(x) (mcintosh(x, tax_lev = tax_lev))
  f8 <- function(x) (menhinick(x, tax_lev = tax_lev))
  f9 <- function(x) (pielou(x, tax_lev = tax_lev, base = base))
  f10 <- function(x) (simpson(x, tax_lev = tax_lev))
  f11 <- function(x) (esimpson(x, tax_lev = tax_lev))
  f12 <- function(x) (fisher(x, tax_lev = tax_lev))

  funs <- list(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12)
  indices <- lapply(funs, function(f) f(x))
  # from list to data.frame
  res <- data.frame(t(matrix(unlist(indices), ncol = length(st.names), byrow = T)))
  rownames(res) <- st.names
  colnames(res) <- c(
    "shannon", "berpar", "brill", "invberpar", "invsimpson",
    "margalef", "mcintosh", "menhinick", "pielou",
    "simpson", "esimpson", "fisher"
  )
  res
}
