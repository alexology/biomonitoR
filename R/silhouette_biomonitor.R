#' @importFrom stats hclust cutree

silhouette_biomonitor <- function(x, nmax = NULL, method = method, clusters = NULL){

  data_clust <- hclust(as.dist(x), method = method)
  dist_ma <- as.matrix(x)
  diag(dist_ma) <- NA



  if(is.null(clusters)){
    res_list <- vector(nmax - 1, mode = "list")
    df_temp <- data.frame(k = numeric(), cluster = character(), mean_dist = numeric())

    for(i in 2:nmax){
      x_ct <- cutree(data_clust, i)
      k <- unique(x_ct)
      for(j in 1:length(k)){
        temp_ma <- dist_ma[x_ct == j , x_ct == j, drop = FALSE]
        temp_ma_not <- dist_ma[!  x_ct == j , x_ct == j, drop = FALSE]
        temp_ma_not <- aggregate(temp_ma_not,by = list(x_ct[x_ct != j]), FUN = mean)[, -1, drop = FALSE]
        res <- data.frame(k = i, cluster = rep(j, ncol(temp_ma)), mean_dist = apply(temp_ma, 2, function(x) sum(na.omit(x))/length(na.omit(x))),
                          min_dist = apply(temp_ma_not, 2, min) )
        df_temp <- rbind(df_temp, res)

      }

    }

    df_temp$res <- ifelse(df_temp[,3] != df_temp[,4], (df_temp[,4] - df_temp[,3]) / pmax(df_temp[,3], df_temp[,4]), 0)
    res_def_k <- aggregate(res ~ cluster + k, df_temp, FUN = mean)
    res_def <- aggregate(res ~ k, df_temp, FUN = function(x) mean((x)) )
    list(x_ct = cutree(data_clust, res_def[which.max(res_def$res), 1]), avg_silh_clusters = res_def_k[res_def_k$k == res_def[which.max(res_def$res), 1] , ], avg_silh = res_def[res_def$k == res_def[which.max(res_def$res), 1] , ]  )

  } else {
    res_list <- vector(nmax - 1, mode = "list")
    df_temp <- data.frame(k = numeric(), cluster = character(), mean_dist = numeric())
    x_ct <- clusters
    k <- unique(x_ct)

    for(i in 2:nmax){
      for(j in 1:length(k)){
        temp_ma <- dist_ma[x_ct == j , x_ct == j, drop = FALSE]
        temp_ma_not <- dist_ma[!  x_ct == j , x_ct == j, drop = FALSE]
        temp_ma_not <- aggregate(temp_ma_not,by = list(x_ct[x_ct != j]), FUN = mean)[, -1, drop = FALSE]
        res <- data.frame(k = i, cluster = rep(j, ncol(temp_ma)), mean_dist = apply(temp_ma, 2, function(x) sum(na.omit(x))/length(na.omit(x))),
                          min_dist = apply(temp_ma_not, 2, min) )
        df_temp <- rbind(df_temp, res)

      }

    }

    df_temp$res <- ifelse(df_temp[,3] != df_temp[,4], (df_temp[,4] - df_temp[,3]) / pmax(df_temp[,3], df_temp[,4]), 0)
    res_def_k <- aggregate(res ~ cluster + k, df_temp, FUN = mean)
    res_def <- aggregate(res ~ k, df_temp, FUN = function(x) mean((x)) )
    list(x_ct = cutree(data_clust, res_def[which.max(res_def$res), 1]), avg_silh_clusters = res_def_k[res_def_k$k == res_def[which.max(res_def$res), 1] , ], avg_silh = res_def[res_def$k == res_def[which.max(res_def$res), 1] , ]  )

  }


}
