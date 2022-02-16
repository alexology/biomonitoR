#' aggregate_taxa
#'
#' @description
#' This function aggregates taxa at various taxonomic levels. Preparing data with `aggregate_taxa()` is necessary for
#' further calculations.
#'
#' @param x Result of `as_biomonitor()`.
#' @param FUN The function to be applied for aggregating at higher taxonomic levels.
#' Must be sum for both abundances and presence-absence data.
#' Default to `sum()`.
#'
#' @keywords aggregate_taxa
#'
#' @importFrom stats aggregate
#'
#' @export
#'
#' @seealso [as_biomonitor]
#'
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#'
#' # example for macrophytes
#' data(oglio)
#'
#' oglio_asb <- as_biomonitor(oglio, group = "mf")
#' oglio_agg <- aggregate_taxa(oglio_asb)
#' richness(oglio_agg, tax_lev = "Species")
#' richness(oglio_agg, tax_lev = "Genus")
#' richness(oglio_agg, tax_lev = "Family")
aggregate_taxa <- function(x, FUN = sum) {

  # check if the object x is of class "biomonitoR"
  if (!inherits(x, "asb")) stop("x is not an object created with as_biomonitor")

  class_x <- class(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  # The following part is to aggregate the dataset provided by the user at different taxonomic levels
  tx <- c("Phylum", "Class", "Subclass", "Order", "Family", "Subfamily", "Tribus", "Genus", "Species", "Subspecies", "Taxa")
  stz <- x[!(names(x) %in% tx)]
  stz_n <- names(stz)

  agg_list <- vector(mode = "list", length = 12)
  names(agg_list) <- c(tx, "Tree")

  for (i in seq_along(tx)) {
    temp.agg <- aggregate(stz, by = list(x[, i]), FUN)
    temp.agg[temp.agg$Group.1 == "", "Group.1"] <- "unassigned"
    names(temp.agg) <- c(tx[i], stz_n)
    # name the first column of the Taxa data.frame "Taxon" for consistency with previous function
    if( identical(tx[i], "Taxa")) {
      names(temp.agg)[1] <- "Taxon"
    }
    agg_list[[i]] <- temp.agg
  }


  agg_list$Tree <- x

  class_x[class_x %in% "asb"] <- "biomonitoR"
  class(agg_list) <- class_x
  agg_list
}
