#' bioco
#'
#' @description
#' bioco function provides site-specific biocontamination index (SBCI), abundance contamination index (ACI) and richness contamination index (RCI) at family rank according to those proposed by Arbaciauskas et al. (2008).
#'
#' @param x Result of function `aggregate_taxa()`.
#' @param alien a vector containing the alien taxa name. Only taxa up to family level will be considered.
#' @param dfref reference database.
#' If a default reference database was used use `mi` for macroinvertebrates or `mf` for macrophytes.
#' Use the same reference database you used for `as_biomonitor()` otherwise.
#' @param digits number of decimal places, default to 2 as in the original work of Arbaciauskas et al. (2008).
#'
#' @details Biocontamination of sampling sites was assessed using a site-specific biocontamination index (SBCI) derived from two metrics: abundance contamination index (ACI) and richness contamination index (RCI) at family/ordinal rank. These indices were calculated as:
#' \deqn{ACI = Na/Nt}
#' where Na and Nt are numbers of specimens of alien taxa and total specimens in a sample, respectively, and
#' \deqn{RCI = Ta/Tt}
#' where Ta is the total number of alien families/orders, and Tt is the total number of identified families/orders. With values of ACI and RCI, the site-specific biocontamination index (SBCI) can then be derived from an easy double entry matrix (Arbaciauskas et al. 2008). Five classes of biocontamination ranging from 0 ("no" contamination) to 4 ("severe" contamination) are defined.
#' A global alien species dataset cannot be provided because alien species may vary among different biogeographical regions and countries. For this reason users need to define and properly upload their own alien reference database based on the studied area they want to study.
#' Applications and examples of these indices are available in Cuk et al. (2019) for Croatian rivers and MacNeil et al. 2010 (The Isle of Man).
#'
#' @keywords aggregate_taxa
#'
#' @references Arbaciauskas, K., Semenchenko, V., Grabowski, M., Leuven, R.S.E.W., Paunovic, M., Son, M.O., Csanyi, B., Gumuliauskaite, S., Konopacka, A., Nehring, S., Van der Velde, G., Vezhnovetz, V., Panov, V.E. (2008). Assessment of biocontamination of benthic macroinvertebrate communities in European inland waterways. Aquatic Invasions 3 (2): 211-230.
#' @references Cuk, R., Milisa, M., Atanackovic, A., Dekic, S., Blazekovic, L., & Zganec, K. (2019). Biocontamination of benthic macroinvertebrate assemblages in Croatian major rivers and effects on ecological quality assessment. Knowledge & Management of Aquatic Ecosystems, (420), 11.
#' @references MacNeil, C., Briffa, M., Leuven, R.S., Gell, F.R. and Selman, R. (2010). An appraisal of a biocontamination assessment method for freshwater macroinvertebrate assemblages; a practical way to measure a significant biological pressure? Hydrobiologia 638:151-159
#'
#' @export
#'
#' @importFrom stats as.formula
#'
#' @seealso [aggregate_taxa]
#'
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' # Toy example:
#' alien <- c("Laccobius", "Setodes bulgaricus", "Caenidae")
#' bioco(data_agr, alien = alien, dfref = "mi")
bioco <- function(x, alien = NULL, dfref = NULL, digits = 2) {

  #  check if the object x is of class "biomonitoR"
  classCheck(x)

  # format names with the first letter as capital letter and the others lowercase

  alien <- sapply(alien, capWords, USE.NAMES = FALSE)

  # if dfref = NULL check if x is of class biomonitor mi, mf or fi, otherwise dfref is needed

  if (inherits(x, "custom") & (identical(dfref, "mi") | identical(dfref, "mf"))) (stop("Please provide the dfref you used for asBiomonitor."))

  if (is.null(dfref)) stop("Please set dfref")
  if (identical(dfref, "mi")) (ref <- mi_ref)
  if (identical(dfref, "mf")) (ref <- mf_ref)
  if (is.data.frame(dfref)) (ref <- dfref)

  # check
  if (is.null(alien)) (stop("Please provide a vector containing the names of alien taxa"))

  alien <- trimws(alien)
  st.names <- names(x[[1]][-1])
  DF <- ref[, 1:10]

  # check if the alien vector contains taxa not present in the reference database
  taxa.vec <- as.character(unique(unlist(DF)))
  if (any(!alien %in% taxa.vec)) {
    absent <- alien[!alien %in% taxa.vec]
    if (length(alien) == length(absent)) {
      absent <- paste(absent, collapse = ", ")
      aci <- rep(0, length(st.names))
      res <- data.frame(aci = aci, rci = aci, sbci = aci)
      rownames(res) <- st.names
      arg.names <- as.list(match.call())[-1]
      arg.names <- arg.names[[which(names(arg.names) == "x")]]
      mes <- paste(absent, "are not part of the taxonomic tree of", as.character(arg.names), sep = " ")
      warning(mes)
      return(res)
    } else {
      absent <- paste(absent, collapse = ", ")
      alien <- alien[alien %in% taxa.vec]
      mes <- paste("The following taxa are absent from the reference database: ", absent)
      warning(mes)
    }
  }

  # Position of taxon in the df data.frame
  taxind <- data.frame(row = numeric(), col = numeric())
  for (i in 1:length(alien)) {
    temp <- which(DF == alien[i], arr.ind = T)
    taxind <- rbind(temp, taxind)
  }

  getAlienAll <- c()
  for (i in 1:nrow(taxind)) {
    a <- taxind[i, 1]
    b <- taxind[i, 2]:10
    temp <- as.character(unlist(DF[a, b]))
    getAlienAll <- c(getAlienAll, temp)
  }

  getAlienAll <- unique(getAlienAll)
  getAlienAll <- getAlienAll[getAlienAll != ""]

  x.taxa <- x[["Tree"]]
  x.taxa <- x.taxa[x.taxa[, "Taxa"] %in% getAlienAll, ]
  if (nrow(x.taxa) == 0) {
    aci <- rep(0, length(st.names))

    if (inherits(x, "bin")) {
      aci <- rep(NA, length(st.names))
    }

    res <- data.frame(aci = aci, rci = aci, sbci = aci)
    rownames(res) <- st.names
    res
  } else {
    x.taxa <- x.taxa[, colnames(x.taxa) %in% c("Family", st.names)]
    x.taxa <- aggregate(as.formula(paste(". ~ ", as.name("Family"))), data = x.taxa, FUN = sum)
    abu.alien <- apply(x.taxa[, -1, drop = FALSE], 2, sum)
    tax.alien <- apply(x.taxa[, -1, drop = FALSE], 2, function(x) sum(x > 0))
    aci <- round(abu.alien / abundance(x, tax_lev = "Taxa", unassigned = TRUE), digits)
    rci <- suppressWarnings(round(tax.alien / richness(x, tax_lev = "Family"), digits))
    cl.lim <- c(1, 0.5, 0.2, 0.1, 0.01, 0)
    cl.lab <- c(0:4)
    cl.abu <- cut(aci, cl.lim, cl.lab, right = TRUE, include.lowest = T)
    cl.tax <- cut(rci, cl.lim, cl.lab, right = TRUE, include.lowest = T)
    cl <- data.frame(as.numeric(as.character(cl.abu)), as.numeric(as.character(cl.tax)))
    sbci <- apply(cl, 1, max)
    if (inherits(x, "bin")) {
      aci <- rep(NA, length(st.names))
      sbci <- rep(NA, length(st.names))
    }
    data.frame(aci = aci, rci = rci, sbci = sbci)
  }
}
