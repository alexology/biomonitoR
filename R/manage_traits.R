#' manage_traits
#'
#' @description
#' A function to select traits based on taxonomic distance.
#'
#' @param x Result of `assign_traits()`.
#' @param method Can be `nearest`, `nearest+`, `nearest-`, `nearest+-` and `neareast-+`.
#'  Please see details for further information.
#' @param traceB If `TRUE` it will return a vector containing taxa excluded from the selection process
#'  because they did not meet the requirments of the selection.
#'
#' @details Method `nearest` selects the traits belonging to the nearest taxa irrispective of their position
#' in the taxonomic tree.
#' Method `nearest+` selects the traits belonging to the nearest taxa that have a taxonomic level equal or finer
#' than the target one. Method `nearest-` does the opposite.
#' Method `nearest+-` selects the traits belonging to the nearest taxa giving priority to taxa having
#' taxonomic level equal or finer than the target one. Method `nearest-+` does the opposite.
#'
#' @export
#' @examples
#' data(macro_ex)
#'
#' data_bio <- as_biomonitor(macro_ex)
#' data_agr <- aggregate_taxa(data_bio)
#' data_ts <- assign_traits(data_agr)
#'
#' # select only the nearest traits
#' data_ts_sub <- manage_traits(data_ts, method = "nearest+-")
manage_traits <- function(x, method = "nearest+-", traceB = FALSE) {
  taxa <- x$Taxa
  x.split <- split(x, taxa)


  if (identical(method, "nearest")) (res <- nearest(x.split))
  if (identical(method, "nearest+")) (res <- nearest_pos(x.split))
  if (identical(method, "nearest-")) (res <- nearest_neg(x.split))
  if (identical(method, "nearest+-")) (res <- nearest_pos_neg(x.split))
  if (identical(method, "nearest-+")) (res <- nearest_neg_pos(x.split))


  res <- do.call("rbind", res)
  rownames(res) <- NULL

  if (!traceB) {
    res
  } else {
    missing.taxa <- taxa[!taxa %in% res$Taxa]
    if (length(missing.taxa) == 0) {
      missing.taxa <- "none"
    }
    list(results = res, taxa_excluded = missing.taxa)
  }
}


nearest <- function(x) {
  temp.res <- lapply(x, FUN = function(z) z[abs(z[, "Taxonomic_distance"]) >= 0, ])
  more_than_0 <- lapply(temp.res, FUN = function(z) (nrow(z[abs(z[, "Taxonomic_distance"]) >= 0, ])))
  if (sum(unlist(more_than_0)) == 0) (stop("Something went wrong. Please contact the mainteiner."))
  temp.res <- temp.res[more_than_0 > 0]
  lapply(temp.res, function(z) z[z[, "Taxonomic_distance"] == min(z[, "Taxonomic_distance"]), ])
}

nearest_pos <- function(x) {
  temp.res <- lapply(x, FUN = function(z) z[z[, "Taxonomic_distance"] >= 0, ])
  more_than_0 <- lapply(temp.res, FUN = function(z) (nrow(z[abs(z[, "Taxonomic_distance"]) >= 0, ])))
  if (sum(unlist(more_than_0)) == 0) (stop("No traits assigned at lower level"))
  temp.res <- temp.res[more_than_0 > 0]
  lapply(temp.res, function(z) z[z[, "Taxonomic_distance"] == min(z[, "Taxonomic_distance"]), ])
}

nearest_neg <- function(x) {
  temp.res <- lapply(x, FUN = function(z) z[z[, "Taxonomic_distance"] <= 0, ])
  more_than_0 <- lapply(temp.res, FUN = function(z) (nrow(z[abs(z[, "Taxonomic_distance"]) >= 0, ])))
  if (sum(unlist(more_than_0)) == 0) (stop("No traits assigned at higher level"))
  temp.res <- temp.res[more_than_0 > 0]
  lapply(temp.res, function(z) z[z[, "Taxonomic_distance"] == max(z[, "Taxonomic_distance"]), ])
}

nearest_pos_neg <- function(x) {
  temp_pos.res <- lapply(x, FUN = function(z) nrow(z[z[, "Taxonomic_distance"] >= 0, ]) > 0)
  temp_neg.res <- lapply(x, FUN = function(z) nrow(z[z[, "Taxonomic_distance"] < 0, ]) > 0)

  if (sum(unlist(temp_pos.res)) == 0 & sum(unlist(temp_neg.res)) == 0) (stop("Something went wrong. Please contact the mainteiner."))


  # assing the nearest traits, at first to positive and then to negative
  if (sum(unlist(temp_pos.res)) > 0) {
    x[unlist(temp_pos.res)] <- lapply(x[unlist(temp_pos.res)], FUN = function(z) z[z[, "Taxonomic_distance"] == min(z[, "Taxonomic_distance"]), ])
  }

  if (sum(unlist(temp_neg.res)) > 0) {
    x[unlist(temp_neg.res)] <- lapply(x[unlist(temp_neg.res)], FUN = function(z) z[z[, "Taxonomic_distance"] == max(z[, "Taxonomic_distance"]), ])
  }

  x
}


nearest_neg_pos <- function(x) {
  temp_neg.res <- lapply(x, FUN = function(z) nrow(z[z[, "Taxonomic_distance"] <= 0, ]) > 0)
  temp_pos.res <- lapply(x, FUN = function(z) nrow(z[z[, "Taxonomic_distance"] > 0, ]) > 0)


  if (sum(unlist(temp_pos.res)) == 0 & sum(unlist(temp_neg.res)) == 0) (stop("Something went wrong. Please contact the mainteiner."))


  # assing the nearest traits, at first to positive and then to negative
  if (sum(unlist(temp_pos.res)) > 0) {
    x[unlist(temp_pos.res)] <- lapply(x[unlist(temp_pos.res)], FUN = function(z) z[z[, "Taxonomic_distance"] == min(z[, "Taxonomic_distance"]), ])
  }

  if (sum(unlist(temp_neg.res)) > 0) {
    x[unlist(temp_neg.res)] <- lapply(x[unlist(temp_neg.res)], FUN = function(z) z[z[, "Taxonomic_distance"] == max(z[, "Taxonomic_distance"]), ])
  }

  x
}
