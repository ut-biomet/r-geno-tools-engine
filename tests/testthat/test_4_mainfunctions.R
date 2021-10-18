# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test of main functions

capture.output({

  # run GWAS ---
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
                               outFile = tempfile(fileext = ".json"))
    },NA)
    expect_error({
      gwas <- readGWAS(gwas_results$file)
    }, NA)
    expect_true(class(gwas) == "list")
    expect_true(all.equal(names(gwas), c("gwas", "metadata")))
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


  # Test run_gwas() with wrong parameters:
  goodParams <- list(genoFile = "../../data/geno/testMarkerData01.vcf.gz",
                     phenoFile = "../../data/pheno/testPhenoData01.csv",
                     genoUrl = NULL,
                     phenoUrl = NULL,
                     trait = "Flowering.time.at.Arkansas",
                     test = "wald",
                     fixed = 0,
                     response = "quantitative",
                     thresh_maf = 0.00,
                     thresh_callrate = 0.99,
                     outFile = NULL)

  wrongParamsL <- list(genoFile = c("/geno-do-not-exist", NA),
                       phenoFile = c("/pheno-do-not-exist", NA),
                       trait = list("Trait.that.do.not.exist", c("a", "b"), NA),
                       test = list("not.a.test", c("c", "d"), NA),
                       fixed = list(-1, 1.5, "1", c(0,1)),
                       response = list("wrong resp", c("e", "f")),
                       thresh_maf = list(-0.06, 0.55, "0.1", c(0, 0.1), NA),
                       thresh_callrate = list(-0.02, 1.42, "0.21", c(0.85, 0.2), NA),
                       outFile = list(c("f1", "f2")))

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("run_gwas, WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          run_gwas(genoFile = params[["genoFile"]],
                   phenoFile = params[["phenoFile"]],
                   trait = params[["trait"]],
                   test = params[["test"]],
                   fixed = params[["fixed"]],
                   response = params[["response"]],
                   thresh_maf = params[["thresh_maf"]],
                   thresh_callrate = params[["thresh_callrate"]],
                   outFile = params[["outFile"]])
        })
      })
    }
  }


    # draw_manhattanPlot ---
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
    expect_error({
      tmpF <- tempfile(fileext = ".html")
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = NA,
                              outFile = tmpF)
    }, NA)
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)

    expect_error({
      tmpF <- tempfile(fileext = ".png")
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = NA,
                              interactive = FALSE,
                              outFile = tmpF)
    }, NA)
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)
  })



  # Test draw_manhattanPlot() with wrong parameters:
  goodParams <- list(gwasFile = "../../data/results/gwasResult.json",
                     gwasUrl = NULL,
                     adj_method = "bonferroni",
                     thresh_p = 0.05,
                     chr = NA,
                     interactive = TRUE,
                     outFile = NULL)

  wrongParamsL <- list(gwasFile = c("doNotExist", NA),
                       gwasUrl = c("doNotExist", NA),
                       adj_method = c("doNotExist", NA),
                       thresh_p = list(-1, 1.02),
                       chr = c("doNotExist", 50),
                       interactive = c('t'),
                       outFile = list(c("f1", "f2")))
  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("draw_manhattanPlot, WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          draw_manhattanPlot(gwasFile = params[["gwasFile"]],
                             gwasUrl = params[["gwasUrl"]],
                             adj_method = params[["adj_method"]],
                             thresh_p = params[["thresh_p"]],
                             chr = params[["chr"]],
                             interactive = params[["interactive"]],
                             outFile = params[["outFile"]])
        })
      })
    }
  }

  # Draw LD plot ----
  test_that("Draw LD plot", {
    expect_error({
      imgFile <- draw_ldPlot(genoFile = "../../data/geno/testMarkerData01.vcf.gz",
                             genoUrl = NULL,
                             from = 42,
                             to = 62,
                             outFile = tempfile(fileext = ".png"))
    },NA)
  })


  # Test draw_ldPlot() with wrong parameters:
  goodParams <- list(genoFile = "../../data/geno/testMarkerData01.vcf.gz",
                     genoUrl = NULL,
                     from = 42,
                     to = 62,
                     outFile = tempfile(fileext = ".png"))

  wrongParamsL <- list(genoFile = c("doNotExist", NA),
                       genoUrl =  c("doNotExist", NA),
                       from = c("42", 42.1),
                       to = c("50", 52.1),
                       outFile = list(c("f1", "f2")))

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("draw_ldPlot, WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          draw_ldPlot(genoFile = params[["genoFile"]],
                      genoUrl = params[["genoUrl"]],
                      from = params[["from"]],
                      to = params[["to"]],
                      outFile = params[["outFile"]])
        })
      })
    }
  }

  # adjust p-values ----
  test_that("GWAS Results adjust p-values", {
    expect_error({
      gwasAdjResults <- run_resAdjustment(gwasFile = "../../data/results/gwasResult.json",
                                          gwasUrl = NULL,
                                          adj_method = "bonferroni",
                                          outFile = tempfile())
    }, NA)
    expect_true(class(gwasAdjResults$gwasAdjusted) == "json")
    expect_error({
      gwas <- readGWAS(gwasAdjResults$file)
    }, NA)
    expect_true(class(gwas) == "list")
    expect_true(all.equal(names(gwas),  c("gwas", "metadata")))
    expect_true(class(gwas$gwas) == "data.frame")
    origGWAS <- readGWAS("../../data/results/gwasResult.json")
    expect_equal(nrow(gwas$gwas), nrow(origGWAS$gwas))
    expect_equal(ncol(gwas$gwas), ncol(origGWAS$gwas)+1)
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

  # Test run_resAdjustment() with wrong parameters:
  goodParams <- list(gwasFile = "../../data/results/gwasResult.json",
                     gwasUrl = NULL,
                     adj_method = "bonferroni",
                     outFile = NULL)

  wrongParamsL <- list(gwasFile = c("doNotExist", NA),
                       gwasUrl = c("doNotExist", NA),
                       adj_method = c("doNotExist", NA),
                       outFile = list(c("f1", "f2")))


  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("run_resAdjustment, WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          run_resAdjustment(gwasFile = params[["gwasFile"]],
                            gwasUrl = params[["gwasUrl"]],
                            adj_method = params[["adj_method"]],
                            outFile = params[["outFile"]])
        })
      })
    }
  }
})
