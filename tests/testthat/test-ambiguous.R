test_that("mcwp_abu", {
  load(ambiguous_df <- system.file("testdata", "ambiguous_df.rda", package="biomonitoR"))
  load(system.file("testdata", "ambiguous_mcwp_abu.rda", package="biomonitoR"))
  load(system.file("testdata", "ambiguous_mcwp_bin.rda", package="biomonitoR"))
  load(system.file("testdata", "ambiguous_rpkc_abu.rda", package="biomonitoR"))
  load(system.file("testdata", "ambiguous_rpkc_bin.rda", package="biomonitoR"))
  ambiguous_asb_abu <- suppressMessages(suppressWarnings(as_biomonitor(ambiguous_df, FUN = sum)))
  ambiguous_agg_abu <- aggregate_taxa(ambiguous_asb_abu)
  ambiguous_asb_bin <- suppressMessages(suppressWarnings(as_biomonitor(ambiguous_df, FUN = bin)))
  ambiguous_agg_bin <- aggregate_taxa(ambiguous_asb_bin)

  expect_equal( suppressWarnings(solve_ambiguous(ambiguous_agg_abu, method = "MCWP")),  ambiguous_mcwp_abu )
  expect_equal( suppressWarnings(solve_ambiguous(ambiguous_agg_bin, method = "MCWP")),  ambiguous_mcwp_bin )
  expect_equal( suppressWarnings(solve_ambiguous(ambiguous_agg_abu, method = "RPKC")),  ambiguous_rpkc_abu )
  expect_equal( suppressWarnings(solve_ambiguous(ambiguous_agg_bin, method = "RPKC")),  ambiguous_rpkc_bin )
})

