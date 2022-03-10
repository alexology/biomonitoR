#' @title transform abundance-coverage data to classes
#'
#' @description From abundance or coverage data it returns classes.
#'
#' @param x A `matrix` or a `data.frame`.
#' @param abundance_classes A numeric vector containing the values for separating classes.
#' @param class_labels A numeric vector containing the labels of the classes. Default to the number of abundance classes plus 1.
#'
#' @details Some indices require abundance or coverage to be transformed in abundance or coverage classes.
#'
#' @export
#'
#' @examples
#'
#' data(mi_prin)
#' abundance_classes(mi_prin)

abundance_classes <- function(x, abundance_classes = c(0.1, 1, 10, 50), class_labels = NULL){

  if(! is.numeric(class_labels) & ! is.null(class_labels)){
    stop("class_labels must be numeric")
  }

  if(is.null(class_labels)){
    class_labels <- 1:(length(abundance_classes) + 1)
  }

  if((length(class_labels) - 1) != length(abundance_classes)){
    stop("length of class_labels must be equal to the length of abundance_classes + 1")
  }


  class.fun <- function(x) cut(x, breaks = c(0, abundance_classes, Inf), labels = class_labels, include.lowest = FALSE, right = TRUE)
  i <- unlist(lapply(x, is.numeric))
  abu.class <- apply(apply(x[, i, drop = FALSE], 2, class.fun), 2, as.numeric)
  abu.class[is.na(abu.class)] <- 0
  x[i] <- abu.class
  x
}
