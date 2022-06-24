test_that("traits", {
  load(system.file("testdata", "macro_traits.rda", package="biomonitoR"))
  load(system.file("testdata", "traits.rda", package="biomonitoR"))
  load(system.file("testdata", "assigned_traits.rda", package="biomonitoR"))

  data_bio <- as_biomonitor(macro_traits, FUN = bin)
  data_agr <- aggregate_taxa(data_bio)
  data_at <- assign_traits(data_agr, trait_db = traits)
  data_at_f <- assign_traits(data_agr, trait_db = traits, filter_by_distance = 0)
  expect_equal(assigned_traits, data_at)
  expect_equal(assigned_traits[assigned_traits$Taxonomic_distance == 0, ], data_at_f)
  expect_equal(assign_traits(data_agr, trait_db = traits, filter_by_distance = "pos"), data_at[data_at[, "Taxonomic_distance"] >= 0, ] )
  expect_equal(assign_traits(data_agr, trait_db = traits, filter_by_distance = "neg"), data_at[data_at[, "Taxonomic_distance"] <= 0, ] )
  expect_error(assign_traits(data_agr, trait_db = traits, filter_by_distance = "eig"), "pos, neg or an integer are needed when filter_by_distance is not NULL" )

  # macrophytes
  # importing data in the biomonitoR format
  oglio_asb <- as_biomonitor(oglio, group = "mf", FUN = bin)
  oglio_agg <- aggregate_taxa(oglio_asb)
  oglio_agg_cust <- oglio_agg
  class(oglio_agg_cust) <- c("biomonitoR", "custom")

  expect_warning(assign_traits(oglio_agg_cust, group = "mf"), "It seems that you used your own reference database. Please check the consistency of the taxonomy used for calculating the index with those of your reference database to have reliable results.")
  expect_error(assign_traits(oglio_agg , trait_db = traits_mf , group = "mf", tax_lev = "Order"), "Maximum taxonomic level is family.")

  # no need to run assign_traits, just to show its use
  oglio_ts <- assign_traits(oglio_agg , trait_db = traits_mf , group = "mf")
  oglio_ts_fam <- assign_traits(oglio_agg , trait_db = traits_mf , group = "mf", tax_lev = "Family")
  expect_equal(dim(oglio_ts), c(44, 28))



})


test_that("sample_traits", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_ts <- assign_traits(data_agr)

  load(system.file("testdata", "sample_traits_run1.rda", package="biomonitoR"))
  load(system.file("testdata", "sample_traits_run2.rda", package="biomonitoR"))

  set.seed(2021)
  expect_equal(sample_traits(data_ts)[, 1, drop = FALSE], sample_traits_run1[, 1, drop = FALSE])
  expect_equal(sample_traits(data_ts)[, 1, drop = FALSE], sample_traits_run2[, 1, drop = FALSE])
})


test_that("fd_indices", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_bio_bin <- suppressWarnings(as_biomonitor(macro_ex, FUN = bin))
  data_agr_bin <- aggregate_taxa(data_bio_bin)
  data_ts <- assign_traits(data_agr)
  data_ts_av <-average_traits(data_ts)
  data_ts_av_m <- data_ts_av
  col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)

  rownames(data_ts_av_m) <- data_ts_av$Taxa
  traits_prep <- ade4::prep.fuzzy(data_ts_av_m[, -1], col.blocks = col_blocks)

  traits_dist <- ade4::ktab.list.df(list(traits_prep))
  traits_dist <- ade4::dist.ktab(traits_dist, type = "F")
  traits_dist_m <- as.matrix(traits_dist)

  taxa_comm <- convert_to_vegan(data_agr, tax_lev = "Taxa")
  taxa_comm <- taxa_comm[, colnames(taxa_comm) %in% rownames(traits_dist_m)]
  traits_dist_m <- traits_dist_m[rownames(traits_dist_m) %in% colnames(taxa_comm), colnames(traits_dist_m) %in% colnames(taxa_comm)]
  trait_dist <- as.dist(traits_dist_m)

  fd_res <- FD::dbFD(trait_dist, taxa_comm, message = FALSE, stand.FRic = FALSE, scale.RaoQ = FALSE, corr = "none")
  check_pcoa_nbdim <- select_pcoa_axes(trait_dist, tresh = 0.91)
  check_pcoa_nbdim_3 <- select_pcoa_axes(trait_dist, tresh = 0.39)

  fd_res_3 <- FD::dbFD(trait_dist, taxa_comm, message = FALSE, stand.FRic = TRUE, scale.RaoQ = TRUE, corr = "none", m = 3)
  f_rich_ex3 <- f_rich(data_agr, trait_db = traits_dist, nbdim = 3)
  f_divs_ex3 <- f_divs(data_agr, trait_db = traits_dist)
  f_red <- f_divs(data_agr, trait_db = traits_dist)
  f_rich_ex3_df <- f_rich(data_agr, trait_db = data_ts_av, type = "F", nbdim = 3, col_blocks = col_blocks)
  f_rich_ex3_df_tb <- f_rich(data_agr, trait_db = data_ts_av, type = "F", nbdim = 3, col_blocks = col_blocks, traceB = TRUE)
  f_eve_ex3 <- f_eve(data_agr, trait_db = traits_dist, nbdim = 33)
  f_red_ex3 <- f_red(data_agr, trait_db = traits_dist)
  f_red_ex3_rao <- f_red_ex3[, 2]
  names(f_red_ex3_rao) <- rownames(f_red_ex3)

  # FD uses all the PCoA dimensions
  fdisp_2 <- FD::fdisp(trait_dist, as.matrix(taxa_comm))[[1]]
  f_disp_ex2 <- f_disp(data_agr, trait_db = traits_dist, nbdim = 33)

  expect_equal(f_rich_ex3, fd_res_3$FRic)
  expect_equal(f_rich_ex3_df, fd_res_3$FRic)
  expect_equal(f_divs_ex3, fd_res_3$RaoQ)
  expect_equal(check_pcoa_nbdim$r2, fd_res$qual.FRic)
  expect_equal(check_pcoa_nbdim_3$r2, fd_res_3$qual.FRic)
  expect_equal(length(f_rich_ex3_df_tb), 7)
  expect_equal(f_red_ex3_rao, fd_res_3$RaoQ)
  expect_equal(f_disp_ex2, fdisp_2)
  expect_equal(f_eve_ex3, fd_res_3$FEve)

  expect_error(f_rich(data_agr), "Please provide trait_db")
  expect_error(f_rich(data_agr, trait_db = "argument"), "trait_db must be a data.frame or a dist object")
  expect_error(f_rich(data_agr, trait_db = data_ts_av, type = NULL), "Please specify a type when trait_db is a data.frame")
  expect_error(f_rich(data_agr, trait_db = data_ts_av, type = 42), "type must be C or F when trait_db is a data.frame")
  expect_error(f_rich(data_agr, trait_db = data_ts_av, type = "C", distance = "gower"), "Using gower distance when type is C is currently not allowed")
  expect_warning(f_rich(data_agr, trait_db = data_ts_av, type = "F", distance = "euclidean", col_blocks = col_blocks), "Are you sure to use euclidean distance when type is F?")
  expect_error(f_rich(data_agr, trait_db = data_ts_av, type = "F"), "Please provide col_blocks")
  expect_error(f_rich(data_agr, trait_db = data_ts_av, type = "F", col_blocks = c(1,1)), "The number of traits in trait_db is not equal to the sum of col_blocks")

  # not really a good example because the trait matrix is euclidean

  f_rich_ex3_cai <- f_rich(data_agr, trait_db = traits_dist, nbdim = 3, correction = "cailliez")
  f_rich_ex3_cai_dist <- suppressWarnings(f_rich(data_agr, trait_db = ade4::cailliez(traits_dist), nbdim = 3, correction = "cailliez"))
  f_rich_ex3_lin <- f_rich(data_agr, trait_db = traits_dist, nbdim = 3, correction = "lingoes")
  f_rich_ex3_lin_dist <- suppressWarnings(f_rich(data_agr, trait_db = ade4::lingoes(traits_dist), nbdim = 3, correction = "cailliez"))
  f_rich_ex3_sqrt <- f_rich(data_agr, trait_db = traits_dist, nbdim = 3, correction = "sqrt")
  f_rich_ex3_sqrt_dist <- suppressWarnings(f_rich(data_agr, trait_db = sqrt(traits_dist), nbdim = 3, correction = "cailliez"))
  f_rich_ex3_quasi <- f_rich(data_agr, trait_db = traits_dist, nbdim = 3, correction = "quasi")
  f_rich_ex3_quasi_dist <- suppressWarnings(f_rich(data_agr, trait_db = ade4::quasieuclid(traits_dist), nbdim = 3, correction = "cailliez"))



  expect_equal(f_rich_ex3_cai_dist, f_rich_ex3_cai_dist)
  expect_equal(f_rich_ex3_lin_dist, f_rich_ex3_lin_dist)
  expect_equal(f_rich_ex3_sqrt_dist, f_rich_ex3_sqrt_dist)
  expect_equal(f_rich_ex3_quasi_dist, f_rich_ex3_quasi_dist)



  data(mi_prin)
  data_bio_non_euclid <- as_biomonitor(mi_prin)
  data_agr_non_euclid <- aggregate_taxa(data_bio_non_euclid)

  data_ts_non_euclid <- assign_traits(data_agr_non_euclid)
  data_ts_av_non_euclid <-average_traits(data_ts_non_euclid)
  data_ts_av_m_non_euclid <- data_ts_av_non_euclid
  col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)

  rownames(data_ts_av_m_non_euclid) <- data_ts_av_non_euclid$Taxa
  traits_prep_non_euclid <- ade4::prep.fuzzy(data_ts_av_m_non_euclid[, -1], col.blocks = col_blocks)

  traits_dist_non_euclid <- ade4::ktab.list.df(list(traits_prep_non_euclid))
  traits_dist_non_euclid <- ade4::dist.ktab(traits_dist_non_euclid, type = "F")
  traits_dist_m_non_euclid <- as.matrix(traits_dist_non_euclid)

  taxa_comm_non_euclid <- convert_to_vegan(data_agr_non_euclid, tax_lev = "Taxa")
  taxa_comm_non_euclid <- taxa_comm_non_euclid[, colnames(taxa_comm_non_euclid) %in% rownames(traits_dist_m_non_euclid)]
  traits_dist_m_non_euclid <- traits_dist_m_non_euclid[rownames(traits_dist_m_non_euclid) %in% colnames(taxa_comm_non_euclid), colnames(traits_dist_m_non_euclid) %in% colnames(taxa_comm_non_euclid)]
  trait_dist_non_euclid <- as.dist(traits_dist_m_non_euclid)



  # 0 distance testing

  trait_dist_0 <- trait_dist_0_test <- as.matrix(trait_dist)
  trait_dist_0[4, 1] <- 0
  trait_dist_0 <- as.dist(trait_dist_0)
  trait_dist_0_test <- as.dist(trait_dist_0_test[rownames(trait_dist_0_test) != "Beraeamyia", colnames(trait_dist_0_test) != "Beraeamyia"])
  comm_0_dist <- macro_ex[macro_ex$Taxa != "Beraeamyia",]
  comm_0_dist[1, 3] <- 4
  data_bio_0_dist <- as_biomonitor(comm_0_dist)
  data_agr_0_dist <- aggregate_taxa(data_bio_0_dist)

  # fake continuous

  traits_continuous <- data_ts_av[, 1:4]
  dist_continuous <- traits_continuous
  rownames(dist_continuous) <- dist_continuous$Taxa
  dist_continuous <- dist(scale(dist_continuous[, -1]))



  ### functional richness -------------------------------------------------------------------------------------


  frich_0_dist <- suppressWarnings(f_rich(data_agr, trait_db = trait_dist_0, nbdim = 3, zerodist_rm = TRUE, traceB = TRUE))
  frich_0_dist_test <- f_rich(data_agr_0_dist, trait_db = trait_dist_0_test, nbdim = 3, traceB = TRUE)


  mess_f_rich <- capture_messages(f_rich(data_agr, trait_db = trait_dist_0, nbdim = 3))
  mess_f_rich <- gsub("\\n", "", mess_f_rich)

  expect_equal(mess_f_rich[1], "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(mess_f_rich[2], "Negative eigenvalues found, please consider to add a correction to the pcoa")


  expect_equal(data.frame(Taxon = "Acentrella", name = "Beraeamyia"), frich_0_dist$duplicated_traits)
  expect_equal(frich_0_dist$taxa, frich_0_dist_test$taxa)
  expect_equal(as.matrix(frich_0_dist$traits), as.matrix(frich_0_dist_test$traits))
  expect_equal(frich_0_dist$results, frich_0_dist_test$results)

  # there is another way to deal with 0, using the option cor.zero of the ade4 functions

  frich_0_dist_cor_zero <- suppressWarnings(f_rich(data_agr, trait_db = trait_dist_0, nbdim = 3, zerodist_rm = FALSE, correction = "cailliez", traceB = TRUE, set_param = list(cor.zero = FALSE)))
  expect_equal(as.matrix(frich_0_dist_cor_zero$traits), suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE))))


  ### continuous traits

  f_rich_cont_df <- suppressMessages(f_rich(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = 2, traceB = TRUE))
  f_rich_cont_ds <- suppressMessages(f_rich(data_agr, trait_db = dist_continuous, nbdim = 2, traceB = TRUE))

  expect_equal(f_rich_cont_df$taxa, f_rich_cont_ds$taxa)
  expect_equal(f_rich_cont_df$results, f_rich_cont_ds$results)
  expect_equal(as.matrix(f_rich_cont_ds$traits), as.matrix(dist_continuous))
  expect_equal(f_rich_cont_df$NA_detection, "No NAs detected")
  expect_equal(f_rich_cont_ds$NA_detection, "NAs cannot be detected when trait_db is a dist object")
  expect_equal(f_rich_cont_df$correction, "none")
  expect_equal(f_rich_cont_ds$correction, "none")
  expect_equal(f_rich_cont_df$nbdim, 2)
  expect_equal(f_rich_cont_ds$nbdim, 2)
  expect_equal(f_rich_cont_df$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(f_rich_cont_ds$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")


  # auto option
  f_rich_cont_df_auto <- suppressMessages(f_rich(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = "auto", traceB = FALSE))
  f_rich_cont_df_3 <- suppressMessages(f_rich(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = 3, traceB = FALSE))

  expect_equal(f_rich_cont_df_auto, f_rich_cont_df_3)



  ### other

    f_rich_dis <- f_rich(data_agr, trait_db = traits_dist, traceB = TRUE)
  f_rich_df <- f_rich(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)


  na_traits <- tidyr::pivot_longer(data_ts_av, -Taxa, names_to = "Traits")
  na_traits <- as.data.frame(na_traits)
  na_traits$Taxa <- as.character(na_traits$Taxa)
  na_traits <- na_traits[is.na(na_traits$value), c("Taxa", "Traits")]
  rownames(na_traits) <- NULL

  f_rich_df_tb <- f_rich(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_rich_dis_bin <- f_rich(data_agr_bin, trait_db = traits_dist)
  f_rich_df_bin <- f_rich(data_agr_bin, trait_db = data_ts_av, type = "F", col_blocks = col_blocks)

  f_rich_dis_cai <- f_rich(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE)
  f_rich_df_cai <- f_rich(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "cailliez", traceB = TRUE)

  f_rich_dis_lin <- f_rich(data_agr_non_euclid, trait_db = lingoes(traits_dist_non_euclid), traceB = TRUE)
  f_rich_df_lin <- f_rich(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "lingoes", traceB = TRUE)

  f_rich_dis_sqrt <- f_rich(data_agr_non_euclid, trait_db = sqrt(traits_dist_non_euclid), traceB = TRUE)
  f_rich_df_sqrt <- f_rich(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "sqrt", traceB = TRUE)

  f_rich_dis_quasi <- f_rich(data_agr_non_euclid, trait_db = quasieuclid(traits_dist_non_euclid), traceB = TRUE)
  f_rich_df_quasi <- f_rich(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "quasi", traceB = TRUE)

  expect_message(f_rich(data_agr_non_euclid, trait_db = traits_dist_non_euclid), "Negative eigenvalues found, please consider to add a correction to the pcoa")
  expect_equal(f_rich_dis$results, f_rich_df$results)
  expect_equal(f_rich_dis_cai$results, f_rich_df_cai$results)
  expect_equal(f_rich_dis_lin$results, f_rich_df_lin$results)
  expect_equal(f_rich_dis_sqrt$results, f_rich_df_sqrt$results)
  expect_equal(f_rich_dis_quasi$results, f_rich_df_quasi$results)
  expect_equal(f_rich_dis_bin, f_rich_df_bin)
  expect_equal(f_rich_df$results, f_rich_df_tb$results)
  expect_error(f_rich(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE, nbdim = "auto", set_param = list(max_nbdim = 9)),
               "there is no optimal number of dimension, please increase the number of dimensions")
  expect_equal(f_rich_dis$duplicated_traits, "no taxa with the same traits")
  expect_equal(f_rich_dis_cai$correction, "none")
  expect_equal(f_rich_dis_lin$correction, "none")
  expect_equal(f_rich_dis_sqrt$correction, "none")
  expect_equal(f_rich_dis_quasi$correction, "none")
  expect_equal(f_rich_df_cai$correction, "cailliez")
  expect_equal(f_rich_df_lin$correction, "lingoes")
  expect_equal(f_rich_df_sqrt$correction, "sqrt")
  expect_equal(f_rich_df_quasi$correction, "quasi")
  expect_equal(f_rich_df$NA_detection, na_traits)



  ### functional evenness -----------------------------------------------------------------------------------

  feve_0_dist <- suppressWarnings(f_eve(data_agr, trait_db = trait_dist_0, nbdim = 3, zerodist_rm = TRUE, traceB = TRUE))
  feve_0_dist_test <- f_eve(data_agr_0_dist, trait_db = trait_dist_0_test, nbdim = 3, traceB = TRUE)


  mess_f_eve <- capture_messages(f_eve(data_agr, trait_db = trait_dist_0, nbdim = 3))
  mess_f_eve <- gsub("\\n", "", mess_f_eve)

  expect_equal(mess_f_eve[1], "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(mess_f_eve[2], "Negative eigenvalues found, please consider to add a correction to the pcoa")


  expect_equal(data.frame(Taxon = "Acentrella", name = "Beraeamyia"), feve_0_dist$duplicated_traits)
  expect_equal(feve_0_dist$taxa, feve_0_dist_test$taxa)
  expect_equal(as.matrix(feve_0_dist$traits), as.matrix(feve_0_dist_test$traits))
  expect_equal(feve_0_dist$results, feve_0_dist_test$results)

  # there is another way to deal with 0, using the option cor.zero of the ade4 functions

  feve_0_dist_cor_zero <- suppressWarnings(f_eve(data_agr, trait_db = trait_dist_0, nbdim = 3, zerodist_rm = FALSE, correction = "cailliez", traceB = TRUE, set_param = list(cor.zero = FALSE)))
  expect_equal(as.matrix(feve_0_dist_cor_zero$traits), suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE))))


  ### continuous traits

  f_eve_cont_df <- suppressMessages(f_eve(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = 2, traceB = TRUE))
  f_eve_cont_ds <- suppressMessages(f_eve(data_agr, trait_db = dist_continuous, nbdim = 2, traceB = TRUE))

  expect_equal(f_eve_cont_df$taxa, f_eve_cont_ds$taxa)
  expect_equal(f_eve_cont_df$results, f_eve_cont_ds$results)
  expect_equal(as.matrix(f_eve_cont_ds$traits), as.matrix(dist_continuous))
  expect_equal(f_eve_cont_df$NA_detection, "No NAs detected")
  expect_equal(f_eve_cont_ds$NA_detection, "NAs cannot be detected when trait_db is a dist object")
  expect_equal(f_eve_cont_df$correction, "none")
  expect_equal(f_eve_cont_ds$correction, "none")
  expect_equal(f_eve_cont_df$nbdim, 2)
  expect_equal(f_eve_cont_ds$nbdim, 2)
  expect_equal(f_eve_cont_df$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(f_eve_cont_ds$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")


  # auto option
  f_eve_cont_df_auto <- suppressMessages(f_eve(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = "auto", traceB = FALSE))
  f_eve_cont_df_3 <- suppressMessages(f_eve(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = 3, traceB = FALSE))

  expect_equal(f_eve_cont_df_auto, f_eve_cont_df_3)



  ### other

  f_eve_dis <- f_eve(data_agr, trait_db = traits_dist, traceB = TRUE)
  f_eve_df <- f_eve(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_eve_df_tb <- f_eve(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_eve_dis_bin <- f_eve(data_agr_bin, trait_db = traits_dist)
  f_eve_df_bin <- f_eve(data_agr_bin, trait_db = data_ts_av, type = "F", col_blocks = col_blocks)

  f_eve_dis_cai <- f_eve(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE)
  f_eve_df_cai <- f_eve(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "cailliez", traceB = TRUE)

  f_eve_dis_lin <- f_eve(data_agr_non_euclid, trait_db = lingoes(traits_dist_non_euclid), traceB = TRUE)
  f_eve_df_lin <- f_eve(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "lingoes", traceB = TRUE)

  f_eve_dis_sqrt <- f_eve(data_agr_non_euclid, trait_db = sqrt(traits_dist_non_euclid), traceB = TRUE)
  f_eve_df_sqrt <- f_eve(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "sqrt", traceB = TRUE)

  f_eve_dis_quasi <- f_eve(data_agr_non_euclid, trait_db = quasieuclid(traits_dist_non_euclid), traceB = TRUE)
  f_eve_df_quasi <- f_eve(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "quasi", traceB = TRUE)

  expect_message(f_eve(data_agr_non_euclid, trait_db = traits_dist_non_euclid), "Negative eigenvalues found, please consider to add a correction to the pcoa")
  expect_equal(f_eve_dis$results, f_eve_df$results)
  expect_equal(f_eve_dis_cai$results, f_eve_df_cai$results)
  expect_equal(f_eve_dis_lin$results, f_eve_df_lin$results)
  expect_equal(f_eve_dis_sqrt$results, f_eve_df_sqrt$results)
  expect_equal(f_eve_dis_quasi$results, f_eve_df_quasi$results)
  expect_equal(f_eve_dis_bin, f_eve_df_bin)
  expect_equal(f_eve_df$results, f_eve_df_tb$results)
  expect_error(f_eve(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE, nbdim = "auto", set_param = list(max_nbdim = 9)),
               "there is no optimal number of dimension, please increase the number of dimensions")
  expect_equal(f_eve_dis$duplicated_traits, "no taxa with the same traits")
  expect_equal(f_eve_dis_cai$correction, "none")
  expect_equal(f_eve_dis_lin$correction, "none")
  expect_equal(f_eve_dis_sqrt$correction, "none")
  expect_equal(f_eve_dis_quasi$correction, "none")
  expect_equal(f_eve_df_cai$correction, "cailliez")
  expect_equal(f_eve_df_lin$correction, "lingoes")
  expect_equal(f_eve_df_sqrt$correction, "sqrt")
  expect_equal(f_eve_df_quasi$correction, "quasi")
  expect_equal(f_eve_df$NA_detection, na_traits)


  expect_error(f_eve(data_agr), "Please provide trait_db")
  expect_error(f_eve(data_agr, trait_db = "argument"), "trait_db must be a data.frame or a dist object")
  expect_error(f_eve(data_agr, trait_db = data_ts_av, type = NULL), "Please specify a type when trait_db is a data.frame")
  expect_error(f_eve(data_agr, trait_db = data_ts_av, type = 42), "type must be C or F when trait_db is a data.frame")
  expect_error(f_eve(data_agr, trait_db = data_ts_av, type = "C", distance = "gower"), "Using gower distance when type is C is currently not allowed")
  expect_warning(f_eve(data_agr, trait_db = data_ts_av, type = "F", distance = "euclidean", col_blocks = col_blocks), "Are you sure to use euclidean distance when type is F?")
  expect_error(f_eve(data_agr, trait_db = data_ts_av, type = "F"), "Please provide col_blocks")
  expect_error(f_eve(data_agr, trait_db = data_ts_av, type = "F", col_blocks = c(1,1)), "The number of traits in trait_db is not equal to the sum of col_blocks")



  ### functional diversity -----------------------------------------------------------------------------------

  fdivs_0_dist <- suppressWarnings(f_divs(data_agr, trait_db = trait_dist_0, zerodist_rm = TRUE, traceB = TRUE))
  fdivs_0_dist_test <- f_divs(data_agr_0_dist, trait_db = trait_dist_0_test, traceB = TRUE)

  expect_error(f_divs(data_agr, trait_db = trait_dist_0), "Non euclidean trait distance. Euclidean property is needed. Please use the correction options
             otherwise consider to remove taxa with the same traits.")



  expect_equal(data.frame(Taxon = "Acentrella", name = "Beraeamyia"), fdivs_0_dist$duplicated_traits)
  expect_equal(fdivs_0_dist$taxa, fdivs_0_dist_test$taxa)
  expect_equal(as.matrix(fdivs_0_dist$traits), as.matrix(fdivs_0_dist_test$traits))
  expect_equal(fdivs_0_dist$results, fdivs_0_dist_test$results)

  # there is another way to deal with 0, using the option cor.zero of the ade4 functions

  fdivs_0_dist_cor_zero <- suppressWarnings(f_divs(data_agr, trait_db = trait_dist_0, zerodist_rm = FALSE, correction = "cailliez", traceB = TRUE, set_param = list(cor.zero = FALSE)))
  temp_divs <- suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE)))
  expect_equal(as.matrix(fdivs_0_dist_cor_zero$traits), suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE))))


  ### continuous traits

  f_divs_cont_df <- suppressMessages(f_divs(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", traceB = TRUE))
  f_divs_cont_ds <- suppressMessages(f_divs(data_agr, trait_db = dist_continuous, traceB = TRUE))

  expect_equal(f_divs_cont_df$taxa, f_divs_cont_ds$taxa)
  expect_equal(f_divs_cont_df$results, f_divs_cont_ds$results)
  expect_equal(as.matrix(f_divs_cont_ds$traits), as.matrix(dist_continuous))
  expect_equal(f_divs_cont_df$NA_detection, "No NAs detected")
  expect_equal(f_divs_cont_ds$NA_detection, "NAs cannot be detected when trait_db is a dist object")
  expect_equal(f_divs_cont_df$correction, "none")
  expect_equal(f_divs_cont_ds$correction, "none")
  expect_equal(f_divs_cont_df$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(f_divs_cont_ds$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")


  ### other

  f_divs_dis <- f_divs(data_agr, trait_db = traits_dist, traceB = TRUE)
  f_divs_df <- f_divs(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_divs_df_tb <- f_divs(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_divs_dis_bin <- f_divs(data_agr_bin, trait_db = traits_dist)
  f_divs_df_bin <- f_divs(data_agr_bin, trait_db = data_ts_av, type = "F", col_blocks = col_blocks)

  f_divs_dis_cai <- f_divs(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE)
  f_divs_df_cai <- f_divs(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "cailliez", traceB = TRUE)

  f_divs_dis_lin <- f_divs(data_agr_non_euclid, trait_db = lingoes(traits_dist_non_euclid), traceB = TRUE)
  f_divs_df_lin <- f_divs(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "lingoes", traceB = TRUE)

  f_divs_dis_sqrt <- f_divs(data_agr_non_euclid, trait_db = sqrt(traits_dist_non_euclid), traceB = TRUE)
  f_divs_df_sqrt <- f_divs(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "sqrt", traceB = TRUE)

  f_divs_dis_quasi <- f_divs(data_agr_non_euclid, trait_db = quasieuclid(traits_dist_non_euclid), traceB = TRUE)
  f_divs_df_quasi <- f_divs(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "quasi", traceB = TRUE)

  expect_error(f_divs(data_agr_non_euclid, trait_db = traits_dist_non_euclid), "Non euclidean trait distance. Euclidean property is needed. Please use the correction options
             otherwise consider to remove taxa with the same traits.")
  expect_equal(f_divs_dis$results, f_divs_df$results)
  expect_equal(f_divs_dis_cai$results, f_divs_df_cai$results)
  expect_equal(f_divs_dis_lin$results, f_divs_df_lin$results)
  expect_equal(f_divs_dis_sqrt$results, f_divs_df_sqrt$results)
  expect_equal(f_divs_dis_quasi$results, f_divs_df_quasi$results)
  expect_equal(f_divs_dis_bin, f_divs_df_bin)
  expect_equal(f_divs_df$results, f_divs_df_tb$results)
  expect_equal(f_divs_dis$duplicated_traits, "no taxa with the same traits")
  expect_equal(f_divs_dis_cai$correction, "none")
  expect_equal(f_divs_dis_lin$correction, "none")
  expect_equal(f_divs_dis_sqrt$correction, "none")
  expect_equal(f_divs_dis_quasi$correction, "none")
  expect_equal(f_divs_df_cai$correction, "cailliez")
  expect_equal(f_divs_df_lin$correction, "lingoes")
  expect_equal(f_divs_df_sqrt$correction, "sqrt")
  expect_equal(f_divs_df_quasi$correction, "quasi")
  expect_equal(f_divs_df$NA_detection, na_traits)


  expect_error(f_divs(data_agr), "Please provide trait_db")
  expect_error(f_divs(data_agr, trait_db = "argument"), "trait_db must be a data.frame or a dist object")
  expect_error(f_divs(data_agr, trait_db = data_ts_av, type = NULL), "Please specify a type when trait_db is a data.frame")
  expect_error(f_divs(data_agr, trait_db = data_ts_av, type = 42), "type must be C or F when trait_db is a data.frame")
  expect_error(f_divs(data_agr, trait_db = data_ts_av, type = "C", distance = "gower"), "Using gower distance when type is C is currently not allowed")
  expect_warning(f_divs(data_agr, trait_db = data_ts_av, type = "F", distance = "euclidean", col_blocks = col_blocks), "Are you sure to use euclidean distance when type is F?")
  expect_error(f_divs(data_agr, trait_db = data_ts_av, type = "F"), "Please provide col_blocks")
  expect_error(f_divs(data_agr, trait_db = data_ts_av, type = "F", col_blocks = c(1,1)), "The number of traits in trait_db is not equal to the sum of col_blocks")


  ### functional dispersion -----------------------------------------------------------------------------------

  fdisp_0_dist <- suppressWarnings(f_disp(data_agr, trait_db = trait_dist_0, nbdim = 3, zerodist_rm = TRUE, traceB = TRUE))
  fdisp_0_dist_test <- f_disp(data_agr_0_dist, trait_db = trait_dist_0_test, nbdim = 3, traceB = TRUE)


  mess_f_disp <- capture_messages(f_disp(data_agr, trait_db = trait_dist_0, nbdim = 3))
  mess_f_disp <- gsub("\\n", "", mess_f_disp)

  expect_equal(mess_f_disp[1], "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(mess_f_disp[2], "Negative eigenvalues found, please consider to add a correction to the pcoa")


  expect_equal(data.frame(Taxon = "Acentrella", name = "Beraeamyia"), fdisp_0_dist$duplicated_traits)
  expect_equal(fdisp_0_dist$taxa, fdisp_0_dist_test$taxa)
  expect_equal(as.matrix(fdisp_0_dist$traits), as.matrix(fdisp_0_dist_test$traits))
  expect_equal(fdisp_0_dist$results, fdisp_0_dist_test$results)

  # there is another way to deal with 0, using the option cor.zero of the ade4 functions

  fdisp_0_dist_cor_zero <- suppressWarnings(f_disp(data_agr, trait_db = trait_dist_0, nbdim = 3, zerodist_rm = FALSE, correction = "cailliez", traceB = TRUE, set_param = list(cor.zero = FALSE)))
  expect_equal(as.matrix(fdisp_0_dist_cor_zero$traits), suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE))))


  ### continuous traits

  f_disp_cont_df <- suppressMessages(f_disp(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = 2, traceB = TRUE))
  f_disp_cont_ds <- suppressMessages(f_disp(data_agr, trait_db = dist_continuous, nbdim = 2, traceB = TRUE))

  expect_equal(f_disp_cont_df$taxa, f_disp_cont_ds$taxa)
  expect_equal(f_disp_cont_df$results, f_disp_cont_ds$results)
  expect_equal(as.matrix(f_disp_cont_ds$traits), as.matrix(dist_continuous))
  expect_equal(f_disp_cont_df$NA_detection, "No NAs detected")
  expect_equal(f_disp_cont_ds$NA_detection, "NAs cannot be detected when trait_db is a dist object")
  expect_equal(f_disp_cont_df$correction, "none")
  expect_equal(f_disp_cont_ds$correction, "none")
  expect_equal(f_disp_cont_df$nbdim, 2)
  expect_equal(f_disp_cont_ds$nbdim, 2)
  expect_equal(f_disp_cont_df$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(f_disp_cont_ds$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")


  # auto option
  f_disp_cont_df_auto <- suppressWarnings(suppressMessages(f_disp(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = "auto", traceB = FALSE)))
  f_disp_cont_df_3 <- suppressWarnings(suppressMessages(f_disp(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", nbdim = 3, traceB = FALSE)))

  expect_equal(f_disp_cont_df_auto, f_disp_cont_df_3)



  ### other

  f_disp_dis <- f_disp(data_agr, trait_db = traits_dist, traceB = TRUE)
  f_disp_df <- f_disp(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_disp_df_tb <- f_disp(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_disp_dis_bin <- f_disp(data_agr_bin, trait_db = traits_dist)
  f_disp_df_bin <- f_disp(data_agr_bin, trait_db = data_ts_av, type = "F", col_blocks = col_blocks)

  f_disp_dis_cai <- f_disp(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE)
  f_disp_df_cai <- f_disp(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "cailliez", traceB = TRUE)

  f_disp_dis_lin <- f_disp(data_agr_non_euclid, trait_db = lingoes(traits_dist_non_euclid), traceB = TRUE)
  f_disp_df_lin <- f_disp(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "lingoes", traceB = TRUE)

  f_disp_dis_sqrt <- f_disp(data_agr_non_euclid, trait_db = sqrt(traits_dist_non_euclid), traceB = TRUE)
  f_disp_df_sqrt <- f_disp(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "sqrt", traceB = TRUE)

  f_disp_dis_quasi <- f_disp(data_agr_non_euclid, trait_db = quasieuclid(traits_dist_non_euclid), traceB = TRUE)
  f_disp_df_quasi <- f_disp(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "quasi", traceB = TRUE)

  expect_message(f_disp(data_agr_non_euclid, trait_db = traits_dist_non_euclid), "Negative eigenvalues found, please consider to add a correction to the pcoa")
  expect_equal(f_disp_dis$results, f_disp_df$results)
  expect_equal(f_disp_dis_cai$results, f_disp_df_cai$results)
  expect_equal(f_disp_dis_lin$results, f_disp_df_lin$results)
  expect_equal(f_disp_dis_sqrt$results, f_disp_df_sqrt$results)
  expect_equal(f_disp_dis_quasi$results, f_disp_df_quasi$results)
  expect_equal(f_disp_dis_bin, f_disp_df_bin)
  expect_equal(f_disp_df$results, f_disp_df_tb$results)
  expect_error(f_disp(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE, nbdim = "auto", set_param = list(max_nbdim = 9)),
               "there is no optimal number of dimension, please increase the number of dimensions")
  expect_equal(f_disp_dis$duplicated_traits, "no taxa with the same traits")
  expect_equal(f_disp_dis_cai$correction, "none")
  expect_equal(f_disp_dis_lin$correction, "none")
  expect_equal(f_disp_dis_sqrt$correction, "none")
  expect_equal(f_disp_dis_quasi$correction, "none")
  expect_equal(f_disp_df_cai$correction, "cailliez")
  expect_equal(f_disp_df_lin$correction, "lingoes")
  expect_equal(f_disp_df_sqrt$correction, "sqrt")
  expect_equal(f_disp_df_quasi$correction, "quasi")
  expect_equal(f_disp_df$NA_detection, na_traits)


  expect_error(f_disp(data_agr), "Please provide trait_db")
  expect_error(f_disp(data_agr, trait_db = "argument"), "trait_db must be a data.frame or a dist object")
  expect_error(f_disp(data_agr, trait_db = data_ts_av, type = NULL), "Please specify a type when trait_db is a data.frame")
  expect_error(f_disp(data_agr, trait_db = data_ts_av, type = 42), "type must be C or F when trait_db is a data.frame")
  expect_error(f_disp(data_agr, trait_db = data_ts_av, type = "C", distance = "gower"), "Using gower distance when type is C is currently not allowed")
  expect_warning(f_disp(data_agr, trait_db = data_ts_av, type = "F", distance = "euclidean", col_blocks = col_blocks), "Are you sure to use euclidean distance when type is F?")
  expect_error(f_disp(data_agr, trait_db = data_ts_av, type = "F"), "Please provide col_blocks")
  expect_error(f_disp(data_agr, trait_db = data_ts_av, type = "F", col_blocks = c(1,1)), "The number of traits in trait_db is not equal to the sum of col_blocks")


  ### functional redundancy -----------------------------------------------------------------------------------

  fred_0_dist <- suppressWarnings(f_red(data_agr, trait_db = trait_dist_0, zerodist_rm = TRUE, traceB = TRUE))
  fred_0_dist_test <- f_red(data_agr_0_dist, trait_db = trait_dist_0_test, traceB = TRUE)

  expect_error(f_red(data_agr, trait_db = trait_dist_0), "Non euclidean trait distance. Euclidean property is needed. Please use the correction options
             otherwise consider to remove taxa with the same traits.")



  expect_equal(data.frame(Taxon = "Acentrella", name = "Beraeamyia"), fred_0_dist$duplicated_traits)
  expect_equal(fred_0_dist$taxa, fred_0_dist_test$taxa)
  expect_equal(as.matrix(fred_0_dist$traits), as.matrix(fred_0_dist_test$traits))
  expect_equal(fred_0_dist$results, fred_0_dist_test$results)

  # there is another way to deal with 0, using the option cor.zero of the ade4 functions

  fred_0_dist_cor_zero <- suppressWarnings(f_red(data_agr, trait_db = trait_dist_0, zerodist_rm = FALSE, correction = "cailliez", traceB = TRUE, set_param = list(cor.zero = FALSE)))
  temp_divs <- suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE)))
  expect_equal(as.matrix(fred_0_dist_cor_zero$traits), suppressWarnings(as.matrix(cailliez(trait_dist_0, cor.zero = FALSE))))


  ### continuous traits

  f_red_cont_df <- suppressMessages(f_red(data_agr, trait_db = traits_continuous, type = "C", distance = "euclidean", traceB = TRUE))
  f_red_cont_ds <- suppressMessages(f_red(data_agr, trait_db = dist_continuous, traceB = TRUE))

  expect_equal(f_red_cont_df$taxa, f_red_cont_ds$taxa)
  expect_equal(f_red_cont_df$results, f_red_cont_ds$results)
  expect_equal(as.matrix(f_red_cont_ds$traits), as.matrix(dist_continuous))
  expect_equal(f_red_cont_df$NA_detection, "No NAs detected")
  expect_equal(f_red_cont_ds$NA_detection, "NAs cannot be detected when trait_db is a dist object")
  expect_equal(f_red_cont_df$correction, "none")
  expect_equal(f_red_cont_ds$correction, "none")
  expect_equal(f_red_cont_df$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")
  expect_equal(f_red_cont_ds$duplicated_traits, "At least a pair of species has the same traits. Depending on your needs, this could be an issue.")


  ### other

  f_red_dis <- f_red(data_agr, trait_db = traits_dist, traceB = TRUE)
  f_red_df <- f_red(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_red_df_tb <- f_red(data_agr, trait_db = data_ts_av, type = "F", col_blocks = col_blocks, traceB = TRUE)

  f_red_dis_bin <- f_red(data_agr_bin, trait_db = traits_dist)
  f_red_df_bin <- f_red(data_agr_bin, trait_db = data_ts_av, type = "F", col_blocks = col_blocks)

  f_red_dis_cai <- f_red(data_agr_non_euclid, trait_db = cailliez(traits_dist_non_euclid), traceB = TRUE)
  f_red_df_cai <- f_red(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "cailliez", traceB = TRUE)

  f_red_dis_lin <- f_red(data_agr_non_euclid, trait_db = lingoes(traits_dist_non_euclid), traceB = TRUE)
  f_red_df_lin <- f_red(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "lingoes", traceB = TRUE)

  f_red_dis_sqrt <- f_red(data_agr_non_euclid, trait_db = sqrt(traits_dist_non_euclid), traceB = TRUE)
  f_red_df_sqrt <- f_red(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "sqrt", traceB = TRUE)

  f_red_dis_quasi <- f_red(data_agr_non_euclid, trait_db = quasieuclid(traits_dist_non_euclid), traceB = TRUE)
  f_red_df_quasi <- f_red(data_agr_non_euclid, trait_db = data_ts_av_non_euclid, type = "F", col_blocks = col_blocks, correction = "quasi", traceB = TRUE)

  expect_error(f_red(data_agr_non_euclid, trait_db = traits_dist_non_euclid), "Non euclidean trait distance. Euclidean property is needed. Please use the correction options
             otherwise consider to remove taxa with the same traits.")
  expect_equal(f_red_dis$results, f_red_df$results)
  expect_equal(f_red_dis_cai$results, f_red_df_cai$results)
  expect_equal(f_red_dis_lin$results, f_red_df_lin$results)
  expect_equal(f_red_dis_sqrt$results, f_red_df_sqrt$results)
  expect_equal(f_red_dis_quasi$results, f_red_df_quasi$results)
  expect_equal(f_red_dis_bin, f_red_df_bin)
  expect_equal(f_red_df$results, f_red_df_tb$results)
  expect_equal(f_red_dis$duplicated_traits, "no taxa with the same traits")
  expect_equal(f_red_dis_cai$correction, "none")
  expect_equal(f_red_dis_lin$correction, "none")
  expect_equal(f_red_dis_sqrt$correction, "none")
  expect_equal(f_red_dis_quasi$correction, "none")
  expect_equal(f_red_df_cai$correction, "cailliez")
  expect_equal(f_red_df_lin$correction, "lingoes")
  expect_equal(f_red_df_sqrt$correction, "sqrt")
  expect_equal(f_red_df_quasi$correction, "quasi")
  expect_equal(f_red_df$NA_detection, na_traits)


  expect_error(f_red(data_agr), "Please provide trait_db")
  expect_error(f_red(data_agr, trait_db = "argument"), "trait_db must be a data.frame or a dist object")
  expect_error(f_red(data_agr, trait_db = data_ts_av, type = NULL), "Please specify a type when trait_db is a data.frame")
  expect_error(f_red(data_agr, trait_db = data_ts_av, type = 42), "type must be C or F when trait_db is a data.frame")
  expect_error(f_red(data_agr, trait_db = data_ts_av, type = "C", distance = "gower"), "Using gower distance when type is C is currently not allowed")
  expect_warning(f_red(data_agr, trait_db = data_ts_av, type = "F", distance = "euclidean", col_blocks = col_blocks), "Are you sure to use euclidean distance when type is F?")
  expect_error(f_red(data_agr, trait_db = data_ts_av, type = "F"), "Please provide col_blocks")
  expect_error(f_red(data_agr, trait_db = data_ts_av, type = "F", col_blocks = c(1,1)), "The number of traits in trait_db is not equal to the sum of col_blocks")




})


test_that("add_bias_to_traits", {

  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  data_ts <- assign_traits(data_agr)
  data_ts_av <- average_traits(data_ts)

  # set col_blocks
  col_blocks <- c(8, 7, 3, 9, 4, 3, 6, 2, 5, 3, 9, 8, 8, 5, 7, 5, 4, 4, 2, 3, 8)

  set.seed(42)
  res <- add_bias_to_traits(data_ts_av, fuzzy = TRUE, col_blocks = col_blocks)[1:5, 1:5]

  data_ts_av_temp <- data_ts_av[, -1]

  set.seed(42)
  res_test <- lapply(data_ts_av_temp, function(x) x + abs(rnorm(length(x), sd = 0.001)))
  res_test <- t(do.call("rbind", res_test))
  res_test <- as.data.frame(ade4::prep.fuzzy(as.data.frame(res_test), col.blocks = col_blocks))[1:5, 1:4]

  expect_equal(res[, -1], res_test)
  expect_error(add_bias_to_traits(data_ts_av, fuzzy = TRUE, col_blocks = NULL), "Please set col_blocks")
  expect_warning(add_bias_to_traits(data_ts_av, fuzzy = FALSE, col_blocks = c(1,2)), "col_blocks will be ignored because fuzzy = FALSE")
})
