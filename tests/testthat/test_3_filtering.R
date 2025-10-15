# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test of filtering functions




capture_output({
  # filterGWAS_pAdj ----
  test_that("filterGWAS_pAdj returns gwas unchanged by default", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA), p_adj = c(0.03, 0.15, 0.3, NA))
    result <- filterGWAS_pAdj(gwas)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_pAdj returns gwas unchanged when filter_pAdj is NULL", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA), p_adj = c(0.03, 0.15, 0.3, NA))
    result <- filterGWAS_pAdj(gwas, filter_pAdj = NULL)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_pAdj returns gwas unchanged when filter_pAdj is 1", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA), p_adj = c(0.03, 0.15, 0.3, NA))
    result <- filterGWAS_pAdj(gwas, filter_pAdj = 1)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_pAdj returns gwas unchanged when gwas is empty", {
    gwas <- data.frame(p = numeric(0), p_adj = numeric(0))
    result <- filterGWAS_pAdj(gwas, filter_pAdj = 0.05)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_pAdj warns and returns gwas when p_adj is NULL", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA))
    expect_warning(
      result <- filterGWAS_pAdj(gwas, filter_pAdj = 0.05),
      "p-values haven't been adjusted"
    )
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_pAdj throws error when filter_pAdj < 0", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA), p_adj = c(0.03, 0.15, 0.3, NA))
    expect_engineError(
      filterGWAS_pAdj(gwas, filter_pAdj = -0.1)
    )
  })

  test_that("filterGWAS_pAdj throws error when filter_pAdj > 1", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA), p_adj = c(0.03, 0.15, 0.3, NA))
    expect_engineError(
      filterGWAS_pAdj(gwas, filter_pAdj = 1.5)
    )
  })

  test_that("filterGWAS_pAdj filters rows correctly based on p_adj threshold", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1), p_adj = c(0.03, 0.15, 0.3))
    result <- filterGWAS_pAdj(gwas, filter_pAdj = 0.1)
    expect_equal(nrow(result), 1)
    expect_equal(result$p_adj, 0.03)
  })

  test_that("filterGWAS_pAdj handles NA values in p_adj correctly", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, 0.2), p_adj = c(0.03, NA, 0.3, 0.4))
    result <- filterGWAS_pAdj(gwas, filter_pAdj = 0.35)
    expect_equal(nrow(result), 2)
    expect_false(any(is.na(result$p_adj)))
  })

  test_that("filterGWAS_pAdj warns when all rows are filtered out", {
    gwas <- data.frame(p = c(0.01, 0.05, 0.1, NA), p_adj = c(0.3, 0.5, 0.7, NA))
    expect_warning(
      result <- filterGWAS_pAdj(gwas, filter_pAdj = 0.1),
      "removed all the points"
    )
    expect_equal(nrow(result), 0)
  })



  # filterGWAS_nPoints ----
  test_that("filterGWAS_nPoints returns gwas unchanged by default", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_nPoints(gwas)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_nPoints returns gwas unchanged when filter_nPoints is NULL", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_nPoints(gwas, filter_nPoints = NULL)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_nPoints returns gwas unchanged when gwas is empty", {
    gwas <- data.frame(id = character(0), p = numeric(0))
    result <- filterGWAS_nPoints(gwas, filter_nPoints = 10)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_nPoints throws error when filter_nPoints < 0", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    expect_engineError(
      filterGWAS_nPoints(gwas, filter_nPoints = -5)
    )
  })

  test_that("filterGWAS_nPoints returns gwas unchanged when filter_nPoints >= nrow", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_nPoints(gwas, filter_nPoints = 5)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_nPoints returns top N points with smallest p-values", {
    gwas <- data.frame(id = letters[1:5], p = c(0.1, 0.01, 0.5, 0.05, 0.001))
    result <- filterGWAS_nPoints(gwas, filter_nPoints = 3)
    expected_results <- gwas[gwas$id %in% c("e", "b", "d"), ]
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })

  test_that("filterGWAS_nPoints handles exact boundary case", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_nPoints(gwas, filter_nPoints = 3)
    expect_equal(nrow(result), 3)
  })

  test_that("filterGWAS_nPoints warns when result is empty (edge case)", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    expect_warning(
      {
        result <- filterGWAS_nPoints(gwas, filter_nPoints = 0)
      },
      "removed all the points"
    )
    expect_equal(nrow(result), 0)
  })



  ## filterGWAS_quant ----
  test_that("filterGWAS_quant returns gwas unchanged by default", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_quant(gwas)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_quant returns gwas unchanged when filter_quant is NULL", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_quant(gwas, filter_quant = NULL)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_quant returns gwas unchanged when filter_quant is 1", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    result <- filterGWAS_quant(gwas, filter_quant = 1)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_quant returns gwas unchanged when gwas is empty", {
    gwas <- data.frame(id = character(0), p = numeric(0))
    result <- filterGWAS_quant(gwas, filter_quant = 0.5)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS_quant throws error when filter_quant < 0", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    expect_engineError(
      filterGWAS_quant(gwas, filter_quant = -0.1)
    )
  })

  test_that("filterGWAS_quant throws error when filter_quant > 1", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    expect_engineError(
      filterGWAS_quant(gwas, filter_quant = 1.5)
    )
  })

  test_that("filterGWAS_quant filters to correct quantile", {
    gwas <- data.frame(
      id = letters[1:10],
      p = c(0.1, 0.01, 0.5, 0.05, 0.001, 0.2, 0.3, 0.4, 0.6, 0.7)
    )
    result <- filterGWAS_quant(gwas, filter_quant = 0.5)
    expect_equal(nrow(result), 5)
    expect_true(all(result$p <= quantile(gwas$p, probs = 0.5)))

    expected_results <- filterGWAS_nPoints(gwas, 5)
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })

  test_that("filterGWAS_quant uses floor for quantile calculation", {
    # ie. handles when the "raw quantile" leads to a non integer number of lines
    # here 3.5
    gwas <- data.frame(
      id = letters[1:10],
      p = c(0.1, 0.01, 0.5, 0.05, 0.001, 0.2, 0.3, 0.4, 0.6, 0.7)
    )
    result <- filterGWAS_quant(gwas, filter_quant = 0.35)
    expect_equal(nrow(result), 3) # floor(10 * 0.35) = 3
    expect_true(all(result$p <= quantile(gwas$p, probs = 0.35)))

    expected_results <- filterGWAS_nPoints(gwas, 3)
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })

  test_that("filterGWAS_quant warns when all rows are filtered out", {
    gwas <- data.frame(id = letters[1:3], p = c(0.01, 0.05, 0.1))
    expect_warning(
      result <- filterGWAS_quant(gwas, filter_quant = 0),
      "removed all the points"
    )
    expect_equal(nrow(result), 0)
  })

  test_that("filterGWAS_quant handles very small quantiles", {
    gwas <- data.frame(id = paste0("id", 1:1000), p = 1:1000 / 1000)
    result <- filterGWAS_quant(gwas, filter_quant = 0.001)
    expect_equal(nrow(result), 1)
    expect_true(all(result$p <= quantile(gwas$p, probs = 0.001)))

    expected_results <- filterGWAS_nPoints(gwas, 1)
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })



  # filterGWAS ----

  test_that("filterGWAS with no filters returns unchanged data", {
    gwas <- data.frame(
      id = letters[1:10],
      p = seq(0.01, 0.1, length.out = 10),
      p_adj = seq(0.02, 0.2, length.out = 10)
    )
    result <- filterGWAS(gwas)
    expect_equal(result, gwas)
  })

  test_that("filterGWAS applies pAdj filter correctly", {
    gwas <- data.frame(
      id = letters[1:10],
      p = seq(0.01, 0.1, length.out = 10),
      p_adj = seq(0.02, 0.2, length.out = 10)
    )

    result <- filterGWAS(gwas, filter_pAdj = 0.1)
    expect_true(nrow(result) < nrow(gwas))
    expect_true(all(result$p_adj <= 0.1))

    expected_results <- filterGWAS_pAdj(gwas, 0.1)
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })

  test_that("filterGWAS applies nPoints filter correctly", {
    gwas <- data.frame(
      id = letters[1:10],
      p = seq(0.01, 0.1, length.out = 10),
      p_adj = seq(0.02, 0.2, length.out = 10)
    )
    result <- filterGWAS(gwas, filter_nPoints = 5)
    expect_equal(nrow(result), 5)
    expect_equal(result$id, letters[1:5])

    expected_results <- filterGWAS_nPoints(gwas, 5)
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })

  test_that("filterGWAS applies quant filter correctly", {
    gwas <- data.frame(
      id = letters[1:10],
      p = seq(0.01, 0.1, length.out = 10),
      p_adj = seq(0.02, 0.2, length.out = 10)
    )

    result <- filterGWAS(gwas, filter_quant = 0.5)
    expect_equal(nrow(result), 5)
    expect_true(all(result$p <= quantile(gwas$p, probs = 0.5)))

    expected_results <- filterGWAS_quant(gwas, 0.5)
    expect_equal(
      result[order(result$id, result$p), ],
      expected_results[order(expected_results$id, expected_results$p), ]
    )
  })

  test_that("filterGWAS chains n_points and quant filters in correct order", {
    gwas <- data.frame(
      id = letters[1:10],
      p = c(seq(0.001, 0.9, length.out = 10))
    )
    gwas$p_adj <- gwas$p

    # Here n_points should have no effects (we keep top 50% (ie. 5) of the p values,
    # then keep best 8 p values)
    expected_results <- filterGWAS(gwas, filter_quant = 0.5)
    result <- filterGWAS(gwas, filter_nPoints = nrow(expected_results) + 3, filter_quant = 0.5)

    expect_equal(result, expected_results)
  })

  test_that("filterGWAS chains pAdj and quant filters in correct order", {
    gwas <- data.frame(
      id = letters[1:10],
      p = c(seq(0.001, 0.9, length.out = 10))
    )
    gwas$p_adj <- gwas$p

    # Here filter_pAdj should have no effects (we keep top 50% of the p values,
    # then keep P values < that the 80% quantile)
    result <- filterGWAS(gwas, filter_pAdj = quantile(gwas$p_adj, 0.8), filter_quant = 0.5)
    expected_results <- filterGWAS(gwas, filter_quant = 0.5)

    expect_equal(result, expected_results)
  })


  test_that("filterGWAS handles empty dataframe", {
    gwas <- data.frame(
      id = character(0),
      p = numeric(0),
      p_adj = numeric(0)
    )

    result <- filterGWAS(
      gwas,
      filter_pAdj = 0.05,
      filter_nPoints = 10,
      filter_quant = 0.5
    )

    expect_equal(result, gwas)
  })
})
