# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test of main function

capture.output({
  test_that("Run GWAS", {
    expect_error({
      gwas_results <- run_gwas(genoFile = "../../data/geno/testMarkerData01.vcf.gz",
                               phenoFile = "../../data/pheno/testPhenoData01.csv",
                               genoUrl = NULL,
                               phenoUrl = NULL,
                               trait = "Flowering.time.at.Arkansas",
                               test = "score",
                               fixed = 0,
                               response = "quantitative",
                               thresh_maf = 0.05,
                               thresh_callrate = 0.95,
                               dir = tempdir())
    },NA)
    expect_error({
      gwas <- readGWAS(gwas_results$file)
    }, NA)
    expect_true(class(gwas) == "list")
    expect_true(all.equal(names(gwas),  c("gwas", "metadata")))
    expect_true(class(gwas_results$gwas) == "json")
    expect_true(class(gwas$gwas) == "data.frame")
    expect_true(class(gwas$metadata) == "list")
    expect_true(
      all.equal(names(gwas$metadata), c("genoFP",
                                        "phenoFP",
                                        "trait",
                                        "test",
                                        "fixed",
                                        "response",
                                        "thresh_maf",
                                        "thresh_callrate",
                                        "date"))
    )

  })



  test_that("Draw Manhattan Plot", {
    expect_error({
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = NA,
                              title = "Example of Manhattan Plot")
    }, "unused argument \\(title")
    expect_error({
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = NA)
    }, NA)
    expect_true(all.equal(class(p), c("plotly", "htmlwidget")))
  })





  test_that("Draw LD plot", {
    expect_error({
      imgFile <- draw_ldPlot(genoFile = "../../data/geno/testMarkerData01.vcf.gz",
                             genoUrl = NULL,
                             from = 42,
                             to = 62,
                             dir = tempdir())
    },NA)
  })



  test_that("GWAS Results adjust p-values", {
    expect_error({
      gwasAdjResults <- run_resAdjustment(gwasFile = "../../data/results/gwasResult.json",
                                          gwasUrl = NULL,
                                          adj_method = "bonferroni",
                                          dir = tempdir())
    }, NA)
    expect_true(class(gwasAdjResults$gwasAdjusted) == "json")
    expect_error({
      gwas <- readGWAS(gwasAdjResults$file)
    }, NA)
    expect_true(class(gwas) == "list")
    expect_true(all.equal(names(gwas),  c("gwas", "metadata")))
    expect_true(class(gwas$gwas) == "data.frame")
    expect_true(class(gwas$metadata) == "list")
    expect_true(
      all.equal(names(gwas$metadata), c("genoFP",
                                        "phenoFP",
                                        "trait",
                                        "test",
                                        "fixed",
                                        "response",
                                        "thresh_maf",
                                        "thresh_callrate",
                                        "date",
                                        "adj_method"))
    )
  })

})
