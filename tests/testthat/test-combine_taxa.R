test_that("combine_taxa", {
  data(macro_ex)
  macro_ex_sub <- macro_ex[2:5, ]

  data_bio <- as_biomonitor(macro_ex_sub)
  data_bio_bin <- suppressWarnings(as_biomonitor(macro_ex_sub, FUN = bin))
  data_agr <- aggregate_taxa(data_bio)
  data_agr_bin <- aggregate_taxa(data_bio_bin)
  res <- combine_taxa(data_agr, tax_lev = "Genus")
  res_bin <- suppressWarnings(combine_taxa(data_agr_bin, tax_lev = "Genus"))

  data_fam_bio <- data_agr[["Genus"]]
  data_fam_bio <- data_fam_bio[! data_fam_bio$Genus %in% "unassigned", ]
  data_fam_bio_bin <- to_bin(data_fam_bio)

  temp_df <- temp_df_bin <- data.frame(Taxa = character(0), Sample_1 = numeric(0), Sample_2 = numeric(0))

  comb_taxa <- combn(3,2)

  for( i in 1:ncol(comb_taxa)){
    data_fam_temp <- data_fam_bio[comb_taxa[, i],]
    data_fam_temp_bin <- data_fam_bio_bin[comb_taxa[, i],]
    temp_df[i, 1] <- paste(data_fam_temp[, 1], collapse = "_")
    temp_df[i, 2:3] <- apply(data_fam_temp[, -1], 2, sum)
    temp_df_bin[i, 1] <- paste(data_fam_temp_bin[, 1], collapse = "_")
    temp_df_bin[i, 2:3] <- apply(data_fam_temp_bin[, -1], 2, sum)
  }

  expect_equal(res, temp_df)
  expect_equal(res_bin, temp_df_bin)
  expect_warning(combine_taxa(data_agr_bin, tax_lev = "Genus"), "combine_taxa with presence-absence data can be meaningless")
})
