# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test of filtering functions




capture_output({
  dta <- readData(
    "../../data/geno/testMarkerData01.vcf.gz",
    "../../data/pheno/testPhenoData01.csv"
  )
  resGwas <- gwas(dta,
    "Flowering.time.at.Arkansas",
    "wald",
    fixed = 0,
    response = "quantitative",
    thresh_maf = 0.05,
    thresh_callrate = 0.95
  )
  adj <- adjustPval(resGwas$p, "bonferroni")
  resGwas$p_adj <- adj$p_adj

  for (f_pAdj in c(1, 0.05)) {

  }
  test_that("filtering function default params", {
    expect_no_error({
      filtdta <- filterGWAS(gwas = resGwas)
    })
    expect_identical(filtdta, resGwas)
  })

  filters <- c(0, runif(5, 0.1, 1), 1)
  for (filter in filters) {
    test_that(paste("filtering function filter_pAdj:", filter), {
      if (filter == 0) {
        expect_warning({
          filtdta <- filterGWAS(
            gwas = resGwas,
            filter_pAdj = filter,
            filter_nPoints = Inf,
            filter_quant = 1
          )
        })
      } else {
        expect_no_error({
          filtdta <- filterGWAS(
            gwas = resGwas,
            filter_pAdj = filter,
            filter_nPoints = Inf,
            filter_quant = 1
          )
        })
      }
      if (filter == 0) {
        expect_true(nrow(filtdta) == 0)
      } else {
        expect_true(min(filtdta$p_adj) < filter)
      }
    })
  }

  filters <- c(0, sample.int(nrow(resGwas), 5), nrow(resGwas))
  for (filter in filters) {
    test_that(paste("filtering function filter_nPoints:", filter), {
      if (filter == 0) {
        expect_warning({
          filtdta <- filterGWAS(
            gwas = resGwas,
            filter_pAdj = 1,
            filter_nPoints = filter,
            filter_quant = 1
          )
        })
      } else {
        expect_no_error({
          filtdta <- filterGWAS(
            gwas = resGwas,
            filter_pAdj = 1,
            filter_nPoints = filter,
            filter_quant = 1
          )
        })
      }
      expect_equal(nrow(filtdta), filter)
    })
  }


  filters <- c(0, runif(5, 0, 1), 1)
  for (filter in filters) {
    test_that(paste("filtering function filter_quant:", filter), {
      if (filter == 0 || floor(nrow(resGwas) * filter) == 0) {
        expect_warning({
          filtdta <- filterGWAS(
            gwas = resGwas,
            filter_pAdj = 1,
            filter_nPoints = Inf,
            filter_quant = filter
          )
        })
      } else {
        expect_no_error({
          filtdta <- filterGWAS(
            gwas = resGwas,
            filter_pAdj = 1,
            filter_nPoints = Inf,
            filter_quant = filter
          )
        })
      }
      expect_equal(nrow(filtdta), floor(nrow(resGwas) * filter))
    })
  }
})
