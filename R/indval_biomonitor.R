indval_biomonitor <-  function(x, clusters, BIN = FALSE){

  clusters <- clusters[match(colnames(x)[-1], names(clusters))]
  x_bin <- to_bin(x)
  rownames(x_bin) <- x_bin[, 1]
  x_bin <- t(x_bin[, -1, drop = FALSE])
  np <- aggregate(x_bin, by = list(clusters), mean)
  cl_names <- np[, "Group.1"]
  np <- np[, -1, drop = FALSE]
  npk <- apply(np, 2, sum)

  if(!BIN){x
    rownames(x) <- x[, 1]
    x <- t(x[, -1, drop = FALSE])
    np_a <- aggregate(x, by = list(clusters), mean)
    np_a <- np_a[, -1, drop = FALSE]
    npk_a <- apply(np_a, 2, sum)
    res <- t(sqrt(sweep(np_a, 2, npk_a, "/")*np))
  } else {
    res <- t(sqrt(sweep(np, 2, npk, "/")*np))
  }


  colnames(res) <- paste("cluster_", cl_names, sep = "")


  res <- data.frame(Taxa = rownames(res), res)
  rownames(res) <- NULL
  res
}
