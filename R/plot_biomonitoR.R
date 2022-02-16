#' @title Plot biomonitoR object
#'
#' @description
#' barplot based on clustering methods
#'
#' @param x An object of class biomonitoR generated with `aggregate_taxa`.
#' @param tax_lev Taxonomic level to be plotted.
#' @param method Clustering method.
#' @param relative Plot absolute or relative abundance.
#' @param unassigned Remove unassigned taxa.
#' @param trans Transformation to be applied to abundance data.
#' @param ... Further arguments to be passed to `plot`.
#'
#' @importFrom plotly plot_ly layout
#' @importFrom ade4 dist.binary
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr pivot_longer %>%
#' @importFrom stats hclust
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#'
#' data(macro_ex)
#' data_bio <- as_biomonitor(mi_prin)
#' data_agr <- aggregate_taxa(data_bio)
#'


plot.biomonitoR <- function(x, tax_lev = "Taxa", method = "ward.D2", relative = TRUE, unassigned = FALSE, trans = NULL, ...){

  if(! inherits(x, "biomonitoR")){
    stop("biomonitoR object needed")
  }

  to_plot <- x[[tax_lev]]
  if(inherits(x, "bin")){
    to_plot <- to_bin(to_plot)
  }


  i <- unlist(lapply(to_plot, is.numeric))

  if(! is.null(trans)){
    to_plot[i] <- trans(to_plot[i])
  }

  # avoid RCMD notes
  Taxa <- NULL

  to_plot_t <- t(to_plot[i])
  names(to_plot)[1] <- "Taxa"



  if(inherits(x, "bin")){
    cluster_dist <- dist.binary(to_plot_t, method = 1)
  } else {
    cluster_dist <- bray_curtis(to_plot_t)
  }

  cluster_data <- hclust(cluster_dist, method = method)

  # to_plot <- to_plot[, c(1, cluster_data$order + 1)]

  if(!unassigned){
    to_plot <- to_plot[! to_plot$Taxa %in% "unassigned",]
  }

  if(relative){
    to_plot[i] <- apply(to_plot[i], 2, function(x) x/sum(x))
  }



  to_plot <-to_plot %>%
    pivot_longer(-Taxa, names_to = "Sites", values_to = "abundance")

  to_plot$Sites <- factor(to_plot$Sites, levels = unique(to_plot$Sites)[cluster_data$order])

  n_taxa <- length(unique(to_plot$Taxa))
  taxa_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_taxa)

  fig <- plot_ly(to_plot,
                x = ~Sites,
                y = ~abundance,
                color = ~Taxa,
                colors = taxa_colors,
                type = "bar") %>%
    layout(barmode = 'stack')


  fig
}
