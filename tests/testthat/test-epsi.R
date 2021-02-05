test_that("epsi", {
  data(macro_ex)
  data_bio <- as_biomonitor(macro_ex)
  data_agr <- aggregate_taxa(data_bio)
  res <- epsi(data_agr)
  res_excel_turley <- c(96.848485, 94.339623)
  names(res_excel_turley) <- colnames(macro_ex)[-1]
  expect_equal(res, res_excel_turley)
})


