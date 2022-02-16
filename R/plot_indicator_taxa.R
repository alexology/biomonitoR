#' @title Plot indicator taxa
#'
#' @description
#' A sankey plot of indicator taxa.
#'
#' @param x Result of `aggregate_taxa()`.
#' @param clusters Sites clusters. If `NULL` an automatic clustering procedure will be used. Otherwise, provide a
#' `data.frame` with 2 columns, the first called `sites` where store sites name while the second called `clusters` where
#' store cluster groups. Please avoid groups with numeric values.
#' @param tax_lev Taxonomic level to be plotted.
#' @param nmax Maximum number of clusters when `clusters` is `NULL`.
#' @param thresh Threshold to identify indicator taxa, ranging from 0 to 100.
#' @param method Clustering method if `clusters` is `NULL`.
#' @param parent Parent taxonomic level. It must be coarser than `tax_lev`.
#' @param opacity Set the opacity of the color link.
#' @param color Set the color link.
#'
#' @details Indicator Taxa Analysis (IndVal; Dufrene and Legendre, 1997) is a method to find indicator species and species assemblages characterizing groups of sites.
#' The function `plot_indicator_taxa` implements the group-equalized IndVal for both presence-absence and abundance data (De Caceres and Legendre, 2009). Please consider
#' that `plot_indicator_taxa` is intended to visualize the main trends present in the data. For having more control on the
#' IndVal calculation consider the `indicspecies` package (De Caceres and Legendre, 2009).
#'
#' @references Dufrene, M., & Legendre, P. (1997). Species assemblages and indicator species: the need for a flexible asymmetrical approach. Ecological monographs, 67(3), 345-366.
#' @references De Caceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566-3574.
#' @export
#' @importFrom plotly plot_ly
#' @importFrom tidyr pivot_longer %>%
#' @importFrom dplyr filter
#' @importFrom stats hclust
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#'
#' data(macro_ex)
#' data_bio <- as_biomonitor(mi_prin)
#' data_agr <- aggregate_taxa(data_bio)
#' plot_indicator_taxa(data_agr)


plot_indicator_taxa <- function(x, clusters = NULL, tax_lev = "Taxa", nmax = 10, thresh = 25, method = "ward.D2", parent = "Order", opacity = 0.8, color = NULL){


  if(! inherits(x, "biomonitoR")){
    stop("biomonitoR object needed")
  }

  n <- sum(unlist(lapply(x[["Tree"]], is.numeric)))
  if(is.null(clusters)){
     if(nmax >= n){
      stop("nmax cannot be greater or equal to the number of samples.")
    }

  } else{
    n_cl <- length(unique(clusters[, "clusters"]))

    if(n_cl >= n){
      stop("the number of clusters cannot be greater or equal to the number of samples.")
    }

  }

  # Avoid RCMD notes
  Taxa <- NULL



  if(is.data.frame(parent)){
    if(sum(!unlist(lapply(parent, is.numeric))) > 1){
        stop("only numerical traits are supported.")
    }

    if(ncol(parent) > 2){
      if(any(rowSums(parent[, -1]) != 1 )){
        stop("when fuzzy data are used the sum of each row needs to be equal to 1.")
      }
    }
  }

  if(! is.data.frame(parent)){
    child_n <- which(colnames(x[["Tree"]]) %in% tax_lev)
    parent_n <- which(colnames(x[["Tree"]]) %in% parent)

    if(parent_n >= child_n){
      stop("parent cannot have a finer taxonomic resolution than tax_lev.")
    }
  }


  to_plot <- x[[tax_lev]]
  if(inherits(x, "bin")){
    to_plot <- to_bin(to_plot)
  }


  i <- unlist(lapply(to_plot, is.numeric))
  to_plot_t <- t(to_plot[i])
  names(to_plot)[1] <- "Taxa"

  if(inherits(x, "bin")){
    cluster_dist <- dist.binary(to_plot_t, method = 1)
  } else {
    cluster_dist <- bray_curtis(to_plot_t)
  }

  x_df <- x[[tax_lev]]
  names(x_df)[1] <- "Taxa"
  x_df <- x_df[! x_df$Taxa %in% "unassigned", ]

  if(! is.data.frame(clusters)){
    res_silh <- silhouette_biomonitor(cluster_dist, nmax = nmax, method = method)
    res_ind <- indval_biomonitor(x_df, clusters = res_silh$x_ct, BIN = inherits(x, "bin"))
  } else {
    clusters <- clusters[match(colnames(x_df)[-1], clusters[, "sites"]), ]
    clusters$x_ct <- as.numeric(as.factor(clusters[, "clusters"]))
    temp_cl <- clusters[, "x_ct"]
    names(temp_cl) <- clusters[, "sites"]
    conv_clusters <- clusters[, c("clusters", "x_ct")]
    conv_clusters <- conv_clusters[! duplicated(conv_clusters), ]
    res_silh <- silhouette_biomonitor(cluster_dist, nmax = nmax, method = method, clusters = temp_cl)
    res_ind <- indval_biomonitor(x_df, clusters = temp_cl, BIN = inherits(x, "bin"))
    cl_res <- names(res_ind)[-1]
    cl_res <- as.numeric(gsub("cluster_", "", cl_res))
    cl_res <- conv_clusters[match(cl_res, conv_clusters[, "x_ct"]),]
    colnames(res_ind)[-1] <- cl_res[, "clusters"]
  }

  clu_meth <- paste("clustering method:", method)
  if(is.null(clusters)){
    clu_num <- paste("best number of clusters:", unique(res_silh$avg_silh$k))
  }



  res_ind <- res_ind %>%
    pivot_longer(-Taxa) %>%
    filter(value >= thresh / 100) %>%
    as.data.frame()

  if(nrow(res_ind) == 0){
    stop("no indicator taxa were found at the selected threshold value.")
  }

  res_ind$value <- round(res_ind$value, 2)
  res_ind$cluster_taxa <- paste(res_ind$name, res_ind$Taxa, sep = " - ")
  res_ind <- res_ind[order(res_ind[, "value"]), ]





  if(! is.data.frame(parent)){
    x_tree <- x[["Tree"]]
    x_tree <- x_tree[, colnames(x_tree) %in% c(parent, tax_lev)]
    x_tree <- x_tree[! duplicated(x_tree),]
    if(any(x_tree[, parent] %in% "")){
      x_tree[x_tree[, parent] %in% "", parent] <- "missing"
    }
    names(x_tree)[2] <- "Taxa"
    res_ind <- merge(res_ind, x_tree, by = "Taxa")
    res_ind$parent_num <- as.numeric(as.factor(res_ind[, parent])) - 1
    res_ind$name_num <- as.numeric(as.factor(res_ind$name)) + max(res_ind$parent_num)
    res_ind$Taxa_num <- as.numeric(as.factor(res_ind$Taxa)) + max(res_ind$name_num)
    res_ind <- res_ind[order(res_ind$parent_num, res_ind$name_num, res_ind$Taxa_num), ]
    lab1 <- unique(res_ind[, parent])
    lab2 <- unique(res_ind$name)
    lab3 <- unique(res_ind$Taxa)
    lab1 <- lab1[order(lab1)]
    lab2 <- lab2[order(lab2)]
    lab3 <- lab3[order(lab3)]

    label <- c(lab1, lab2, lab3)
    source <- c(res_ind$parent_num, res_ind$name_num)
    target <- c(res_ind$name_num, res_ind$Taxa_num)
    value <- c(res_ind$value, res_ind$value ) * 100

    n_clus <- length(lab1)
    clus_colors <- colorRampPalette(brewer.pal(8, "Blues"))(n_clus)
    n_taxa <- length(lab2)
    taxa_colors <- colorRampPalette(brewer.pal(8, "Oranges"))(n_taxa)
    n_parent <- length(lab3)
    parent_colors <- colorRampPalette(brewer.pal(8, "Greens"))(n_parent)

    res_ind_lab1 <- res_ind[ , c(parent, "value")]
    colnames(res_ind_lab1)[1] <- "traits"
    res_ind_lab1 <- aggregate(value ~ traits, res_ind_lab1, FUN = mean)
    res_ind_lab1 <- res_ind_lab1[order(res_ind_lab1$traits),]
    res_ind_lab1 <- round(res_ind_lab1[, "value"] * 100, 0)
    res_ind_lab2 <- res_ind[ , c("name", "value")]
    res_ind_lab2 <- aggregate(value ~ name, res_ind_lab2, FUN = mean)
    res_ind_lab2 <- res_ind_lab2[order(res_ind_lab2$name),]
    res_ind_lab2 <- round(res_ind_lab2[, "value"] * 100, 0)
    res_ind_lab3 <- res_ind[ , c("Taxa", "value")]
    res_ind_lab3 <- res_ind_lab3[! duplicated(res_ind_lab3), ]
    res_ind_lab3 <- res_ind_lab3[order(res_ind_lab3$Taxa),]
    res_ind_lab3 <- res_ind_lab3[, "value"] * 100
    txt <- c(res_ind_lab1, res_ind_lab2, res_ind_lab3)



  } else {
    x_tree <- parent %>%
      pivot_longer(-Taxa, names_to = "traits", values_to = "traits_value")

    res_ind <- merge(res_ind, x_tree, by = "Taxa")
    res_ind <- res_ind[res_ind$traits_value > 0 ,]
    res_ind$value_traits_iv <- res_ind$value * res_ind$traits_value
    res_ind$traits_num <- as.numeric(as.factor(res_ind$traits)) - 1
    res_ind$name_num <- as.numeric(as.factor(res_ind$name)) + max(res_ind$traits_num)
    res_ind$Taxa_num <- as.numeric(as.factor(res_ind$Taxa)) + max(res_ind$name_num)
    res_ind <- res_ind[order(res_ind$traits_num, res_ind$name_num, res_ind$Taxa_num), ]
    res_ind_sub <- res_ind[, c("name", "Taxa", "name_num", "Taxa_num", "value")]
    res_ind_sub <- res_ind_sub[!duplicated(res_ind_sub), ]
    lab1 <- unique(res_ind$traits)
    lab2 <- unique(res_ind$name)
    lab3 <- unique(res_ind$Taxa)
    lab1 <- lab1[order(lab1)]
    lab2 <- lab2[order(lab2)]
    lab3 <- lab3[order(lab3)]

    label <- c(lab1, lab2, lab3)
    source <- c(res_ind$traits_num, res_ind_sub$name_num)
    target <- c(res_ind$name_num, res_ind_sub$Taxa_num)
    value <- c(res_ind$value_traits_iv, res_ind_sub$value ) * 100

    n_clus <- length(lab1)
    clus_colors <- colorRampPalette(brewer.pal(8, "Blues"))(n_clus)
    n_taxa <- length(lab2)
    taxa_colors <- colorRampPalette(brewer.pal(8, "Oranges"))(n_taxa)
    n_parent <- length(lab3)
    parent_colors <- colorRampPalette(brewer.pal(8, "Greens"))(n_parent)

    res_ind_lab1 <- res_ind[ , c("traits", "value")]
    res_ind_lab1 <- aggregate(value ~ traits, res_ind_lab1, FUN = mean)
    res_ind_lab1 <- round(res_ind_lab1[, "value"] * 100, 0)
    res_ind_lab2 <- res_ind[ , c("name", "value")]
    res_ind_lab2 <- aggregate(value ~ name, res_ind_lab2, FUN = mean)
    res_ind_lab2 <- round(res_ind_lab2[, "value"] * 100, 0)
    res_ind_lab3 <- res_ind[ , c("Taxa", "value")]
    res_ind_lab3 <- res_ind_lab3[! duplicated(res_ind_lab3), "value"] * 100
    txt <- c(res_ind_lab1, res_ind_lab2, res_ind_lab3)
  }


  fig <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = " IV",


    node = list(
      label = label,
      customdata = txt,
      color = c(clus_colors, taxa_colors, parent_colors),
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      ),
      hovertemplate = "%{label}<br>IndVal: %{customdata}<extra></extra>"
    ),

    link = list(
      source = source,
      target = target,
      value =  value,
      color = paste("rgba(211, 211, 211, ", opacity, ")")

    )
  )

  fig <- fig %>% layout(
    title = "indicator taxa analysis",
    font = list(
      size = 10
    ),
    xaxis = list(showgrid = FALSE, zeroline = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE)
  )

  if(is.null(clusters)){
    clu_iv_num <- paste("number of indicator taxa:", length(unique(res_ind$Taxa)))
    clu_num_def <- paste("number of clusters with indicator taxa:",  length(unique(res_ind$name)))
    clu_avg_silh <-  paste("average silhouette:", round(res_silh$avg_silh[2], 2))
    cat(paste(clu_meth, clu_num, clu_iv_num, clu_num_def, clu_avg_silh, sep = "\n"))
  }

  fig

}
