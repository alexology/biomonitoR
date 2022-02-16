#' show_scores
#'
#' This function print the databases used for the calculation of several indices implemented in `biomonitoR`.
#' @param index Name of the index for which information are needed.
#' @param method Method of the index specified in `index`.
#' @keywords as_biomonitor
#' @export
#' @seealso \code{\link{aspt}}, \code{\link{bmwp}}, \code{\link{life}}, \code{\link{whpt}}, \code{\link{psi}}, \code{\link{epsi}}
#' @examples
#' show_scores(index = "aspt", method = "spa")
#' show_scores(index = "life", method = "extence")
show_scores <- function(index = NULL, method = NULL) {
  to.print <- file_list[file_list$Index == index & file_list$Method == method, , drop = FALSE]

  if (nrow(to.print) == 0) {
    mes0 <- paste("Method", method, "is not implemented into the index", index, "or index and methods names are wrong", sep = " ")
    mes1 <- paste("Index must be one of: ", paste(as.character(unique(file_list[, "Index"])), collapse = ", "))
    mes2 <- paste("Method must be one of: ", paste(as.character(unique(file_list[, "Method"])), collapse = ", "))
    stop(cat(mes0, mes1, mes2, sep = "\n"))
  }

  to.print.list <- list()

  for (i in 1:nrow(to.print)) {
    to.print.list[[i]] <- get(as.character(to.print[i, "File"]))
  }

  names(to.print.list) <- to.print[, "Type"]
  to.print.list
}
