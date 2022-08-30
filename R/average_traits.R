#' @describeIn assign_traits average traits values for fuzzy data
#' @importFrom ade4 prep.fuzzy.var

average_traits <- function(x, col_blocks = NULL, type = "F") {

  if(! type %in% c("F", "C")){
    stop("type must be one of F or C")
  }

  if(identical(type, "F")){
    if (is.null(col_blocks)) {
      col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)
    } else {
      col_blocks <- col_blocks
    }

    x <- x[, -c(2:5), drop = FALSE]
    x[is.na(x)] <- 0
    trait_mean <- aggregate(. ~ Taxa, data = x, FUN = mean)
    rownames(trait_mean) <- trait_mean[, 1]
    trait_mean <- trait_mean[, -1, drop = FALSE]
    trait_res <- prep.fuzzy(trait_mean, col.blocks = col_blocks)
    rownames(trait_res) <- NULL
    res <- data.frame(Taxa = as.factor(rownames(trait_mean)), trait_res)
  }

  if(identical(type, "C")){
    x <- x[, -c(2:5), drop = FALSE]
    x[is.na(x)] <- 0
    res <- aggregate(. ~ Taxa, data = x, FUN = mean)
    rownames(res) <- NULL
    res[, "Taxa"] <- as.factor(res[, "Taxa"])
  }

  res
}
