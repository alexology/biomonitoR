#' @title Plot objects generated with as_biomonitor
#'
#' @description
#' Plot 2 taxonomic levels using a sunburst visualization.
#'
#' @param x Result of `as_biomonitor()`.
#' @param parent Name of the coarser taxonomic level.
#' @param child Name of the finer taxonomic level.
#' @param type Plot present-absence data with `pa`, abundance data with `abundance` and
#' detection probabilities with `frequency`.
#' @param remove_empty_child Remove empty level of the child taxonomic level.
#' @param trans Transformation for abundance data.
#' @param ... Further arguments to be passed to `plot`.
#'
#' @export
#' @importFrom plotly plot_ly
#'
#' @examples
#' data(macro_ex)
#' data_bio <- as_biomonitor(macro_ex)
#' plot(data_bio)

plot.asb <- function(x, parent = "Order", child = "Family", type = "pa", remove_empty_child = FALSE, trans = NULL, ...){

  if(! inherits(x, "asb")){
    stop("asb object needed")
  }

  to_plot <- x[[1]]
  i <- unlist(lapply(to_plot, is.numeric))



  # check if the taxonomic level of parent is greater than those of child
  tax_colnames <- colnames(to_plot[!i])

  if(which(tax_colnames == parent) >= which(tax_colnames == child)){
    stop("The taxonomic level of parent needs to be lower than child.")
  }


  i <- (colnames(to_plot) %in% c(parent, child)) | unlist(lapply(to_plot, is.numeric))
  to_plot <- to_plot[i]
  names(to_plot)[1:2] <- c("Parent", "Child")
  to_plot[to_plot[, "Parent"] == "", "Parent"] <- "missing"

  if(remove_empty_child){
    to_plot<- to_plot[! to_plot[, "Child"] %in% "", ]
  }

  to_plot <- aggregate(. ~ Parent + Child, data = to_plot, FUN = sum)
  to_plot$merged <- paste(to_plot[, "Parent"], to_plot[, "Child"], sep = " - ")

  parent_lab <- unique(to_plot[, "Parent"])

  lab_id <- c(parent_lab, to_plot$merged)
  lab_lab <- c(parent_lab, to_plot[, "Child"])
  lab_lab <- gsub(" - ", "<br>", lab_lab)
  lab_par <- c(rep("", length(parent_lab)), to_plot[, "Parent"])


  if(!identical(type, "pa") & !inherits(x, "bin")){

    if(identical(type, "abundance")){
      i <- unlist(lapply(to_plot, is.numeric))
      to_plot <- data.frame(to_plot[,! i], values = apply(to_plot[, i], 1, sum))
      if(! is.null(trans)){
        to_plot[, "values"] <- do.call(trans, list(to_plot[, "values"]))
      }
    }
    if(identical(type, "frequency")){
      i <- unlist(lapply(to_plot, is.numeric))
      to_plot <- data.frame(to_plot[,! i], values = apply(to_plot[, i], 1, function(x) sum(x > 0)) / ncol(to_plot[, i]) *100)
    }

    agg_res <- aggregate(to_plot[, "values", drop = FALSE], by = list(factor(to_plot[, "Parent"])), FUN = sum)
    agg_res <- agg_res[match(parent_lab, agg_res[, "Group.1"]),]
    lab_val <- c(agg_res[, "values"], to_plot[, "values"])

    fig <- plot_ly(
      ids = lab_id,
      labels = lab_lab,
      parents = lab_par,
      values = lab_val,
      type = 'sunburst',
      branchvalues = 'total')
  }

  if(identical(type, "pa") | inherits(x, "bin")){

    fig <- plot_ly(
      ids = lab_id,
      labels = lab_lab,
      parents = lab_par,
      type = 'sunburst')
  }

  fig

}
