test_that("select_pcoa", {

  data(mi_prin)
  data_bio <- as_biomonitor(mi_prin)
  data_agr <- aggregate_taxa(data_bio)
  data_ts <- assign_traits(data_agr)
  data_ts_av <-average_traits(data_ts)
  col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)

  rownames(data_ts_av) <- data_ts_av$Taxa
  traits_prep <- ade4::prep.fuzzy(data_ts_av[, -1], col.blocks = col_blocks)

  traits_dist <- ade4::ktab.list.df(list(traits_prep))
  traits_dist <- ade4::dist.ktab(traits_dist, type = "F")
  traits_dist_m <- as.matrix(traits_dist)

  taxa_comm <- convert_to_vegan(data_agr, tax_lev = "Taxa")
  taxa_comm <- taxa_comm[, colnames(taxa_comm) %in% rownames(traits_dist_m)]
  traits_dist_m <- traits_dist_m[rownames(traits_dist_m) %in% colnames(taxa_comm), colnames(traits_dist_m) %in% colnames(taxa_comm)]
  trait_dist <- as.dist(traits_dist_m)

  trait_spcoa <- select_pcoa_axes(traits_dist, tresh = 0.2)
  trait_spcoa_sqrt <- select_pcoa_axes(traits_dist, tresh = 0.1)

  fd_res_none <- FD::dbFD(traits_dist, taxa_comm, message = FALSE, stand.FRic = TRUE,  corr = "none", m = 2)
  fd_res_cai <- FD::dbFD(traits_dist, taxa_comm, message = FALSE, stand.FRic = TRUE,  corr = "cailliez", m = 2)
  fd_res_lin <- FD::dbFD(traits_dist, taxa_comm, message = FALSE, stand.FRic = TRUE,  corr = "lingoes", m = 2)
  fd_res_sqrt <- FD::dbFD(traits_dist, taxa_comm, message = FALSE, stand.FRic = TRUE,  corr = "sqrt", m = 2)

  expect_equal(trait_spcoa[1, 4], fd_res_none$qual.FRic)
  expect_equal(trait_spcoa[2, 4], fd_res_cai$qual.FRic)
  expect_equal(trait_spcoa[3, 4], fd_res_lin$qual.FRic)
  expect_equal(trait_spcoa_sqrt[4, 4], fd_res_sqrt$qual.FRic)


  trait_spcoa_cor <- select_pcoa_axes(traits_dist, method = "cor", tresh = 0.7)
  pcoa_none <- suppressWarnings(cmdscale(trait_dist, k = trait_spcoa_cor[1, 3]))
  pcoa_cai <- suppressWarnings(cmdscale(ade4::cailliez(trait_dist), k = trait_spcoa_cor[2, 3]))
  pcoa_lin <- suppressWarnings(cmdscale(ade4::lingoes(trait_dist), k = trait_spcoa_cor[3, 3]))
  pcoa_sqrt <- suppressWarnings(cmdscale(sqrt(trait_dist), k = trait_spcoa_cor[4, 3]))
  pcoa_quasi <- suppressWarnings(cmdscale(ade4::quasieuclid(trait_dist), k = trait_spcoa_cor[5, 3]))

  expect_equal(trait_spcoa_cor[1, 4], cor(dist(pcoa_none), traits_dist))
  expect_equal(trait_spcoa_cor[2, 4], cor(dist(pcoa_cai), traits_dist))
  expect_equal(trait_spcoa_cor[3, 4], cor(dist(pcoa_lin), traits_dist))
  expect_equal(trait_spcoa_cor[4, 4], cor(dist(pcoa_sqrt), traits_dist))
  expect_equal(trait_spcoa_cor[5, 4], cor(dist(pcoa_quasi), traits_dist))


  trait_spcoa_mai <- select_pcoa_axes(traits_dist, method = "maire", tresh = 0.05)

  pcoa_none <- suppressWarnings(cmdscale(trait_dist, k = trait_spcoa_mai[1, 3]))
  pcoa_cai <- suppressWarnings(cmdscale(ade4::cailliez(trait_dist), k = trait_spcoa_mai[2, 3]))
  pcoa_lin <- suppressWarnings(cmdscale(ade4::lingoes(trait_dist), k = trait_spcoa_mai[3, 3]))
  pcoa_sqrt <- suppressWarnings(cmdscale(sqrt(trait_dist), k = trait_spcoa_mai[4, 3]))
  pcoa_quasi <- suppressWarnings(cmdscale(ade4::quasieuclid(trait_dist), k = trait_spcoa_mai[5, 3]))

  S <- nrow(as.matrix(trait_dist))

  dist_none_mai <- dist(pcoa_none, method = "euclidean")
  y_none_mai <- dist_none_mai/ max(dist_none_mai) * max(trait_dist)
  sd_none_mai <- round(((sum((trait_dist - y_none_mai)^2)) / (S * (S - 1) / 2)), 6)

  dist_cai_mai <- dist(pcoa_cai, method = "euclidean")
  y_cai_mai <- dist_cai_mai/ max(dist_cai_mai) * max(ade4::cailliez(trait_dist))
  sd_cai_mai <- round(((sum((ade4::cailliez(trait_dist) - y_cai_mai)^2)) / (S * (S - 1) / 2)), 6)

  dist_lin_mai <- dist(pcoa_lin, method = "euclidean")
  y_lin_mai <- dist_lin_mai/ max(dist_lin_mai) * max(ade4::lingoes(trait_dist))
  sd_lin_mai <- round(((sum((ade4::lingoes(trait_dist) - y_lin_mai)^2)) / (S * (S - 1) / 2)), 6)

  dist_sqrt_mai <- dist(pcoa_sqrt, method = "euclidean")
  y_sqrt_mai <- dist_sqrt_mai/ max(dist_sqrt_mai) * max(sqrt(trait_dist))
  sd_sqrt_mai <- round(((sum((sqrt(trait_dist) - y_sqrt_mai)^2)) / (S * (S - 1) / 2)), 6)

  dist_quasi_mai <- dist(pcoa_quasi, method = "euclidean")
  y_quasi_mai <- dist_quasi_mai/ max(dist_quasi_mai) * max(ade4::quasieuclid(trait_dist))
  sd_quasi_mai <- round(((sum((ade4::quasieuclid(trait_dist) - y_quasi_mai)^2)) / (S * (S - 1) / 2)), 6)

  expect_equal(trait_spcoa_mai[1, 4], sd_none_mai)
  expect_equal(trait_spcoa_mai[2, 4], sd_cai_mai)
  expect_equal(trait_spcoa_mai[3, 4], sd_lin_mai)
  expect_equal(trait_spcoa_mai[4, 4], sd_sqrt_mai)
  expect_equal(trait_spcoa_mai[5, 4], sd_quasi_mai)

})
