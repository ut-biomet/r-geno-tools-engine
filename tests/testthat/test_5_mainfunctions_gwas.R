capture.output({
  # run GWAS ----
  files <- list(
    test_01 = list(
      geno = "../../data/geno/testMarkerData01.vcf.gz",
      pheno = "../../data/pheno/testPhenoData01.csv",
      trait = "Flowering.time.at.Arkansas"
    ),
    test_singleTrait = list(
      geno = "../data/Result_genos_hd_subset-initialColl.vcf.gz",
      pheno = "../data/resistance_initColl.csv",
      trait = "resist"
    ),
    missingSNPid = list(
      geno = "../data/Result_genos_hd_subset-initialColl_missingId.vcf.gz",
      pheno = "../data/resistance_initColl.csv",
      trait = "resist"
    ),
    unused_trait_badly_formatted = list(
      geno = "../data/geno_G1_duplicated_ids.vcf.gz",
      pheno = "../../data/genomic_selection/pheno_train_with_pheno2_trait_badly_formated.csv",
      trait = "pheno"
    )
  )
  for (test in names(files)) {
    test_that(paste0("Run GWAS", test), {
      expect_no_error({
        gwas_results <- run_gwas(
          genoFile = files[[test]]$geno,
          phenoFile = files[[test]]$pheno,
          genoUrl = NULL,
          phenoUrl = NULL,
          trait = files[[test]]$trait,
          test = "score",
          fixed = 0,
          response = "quantitative",
          thresh_maf = 0.05,
          thresh_callrate = 0.95,
          outFile = tempfile(fileext = ".json")
        )
      })
      expect_no_error({
        gwas <- readGWAS(gwas_results$file)
      })
      expect_true(class(gwas) == "list")
      expect_true(all.equal(names(gwas), c("gwas", "metadata")))
      expect_true(class(gwas_results$gwas) == "json")
      expect_true(class(gwas$gwas) == "data.frame")
      expect_true(class(gwas$metadata) == "list")
      expect_true(
        all.equal(names(gwas$metadata), c(
          "genoFP",
          "phenoFP",
          "trait",
          "test",
          "fixed",
          "response",
          "thresh_maf",
          "thresh_callrate",
          "date"
        ))
      )
    })
  }


  # Test run_gwas() with wrong parameters:
  goodParams <- list(
    genoFile = "../../data/geno/testMarkerData01.vcf.gz",
    phenoFile = "../../data/pheno/testPhenoData01.csv",
    genoUrl = NULL,
    phenoUrl = NULL,
    trait = "Flowering.time.at.Arkansas",
    test = "wald",
    fixed = 0,
    response = "quantitative",
    thresh_maf = 0.00,
    thresh_callrate = 0.99,
    outFile = NULL
  )

  wrongParamsL <- list(
    genoFile = c("/geno-do-not-exist", NA),
    phenoFile = c(
      "/pheno-do-not-exist", NA,
      "../data/pheno_duplicated.csv",
      "../data/pheno_no_common_ind_with_geno.csv",
      "../data/pheno_inconsistent_both_numeric_and_character_values.csv"
    ),
    trait = list("Trait.that.do.not.exist", c("a", "b"), NA),
    test = list("not.a.test", c("c", "d"), NA),
    fixed = list(-1, 1.5, "1", c(0, 1)),
    response = list("wrong resp", c("e", "f")),
    thresh_maf = list(-0.06, 0.55, "0.1", c(0, 0.1), NA),
    thresh_callrate = list(-0.02, 1.42, "0.21", c(0.85, 0.2), NA),
    outFile = list(c("f1", "f2"))
  )

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i + 1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("run_gwas, WrongParams", p, i, sep = "-")

      test_that(testName, {
        err <- expect_engineError({
          run_gwas(
            genoFile = params[["genoFile"]],
            phenoFile = params[["phenoFile"]],
            trait = params[["trait"]],
            test = params[["test"]],
            fixed = params[["fixed"]],
            response = params[["response"]],
            thresh_maf = params[["thresh_maf"]],
            thresh_callrate = params[["thresh_callrate"]],
            outFile = params[["outFile"]]
          )
        })
      })
    }
  }


  # draw_manhattanPlot ----
  test_that("Draw Manhattan Plot", {
    expect_no_error({
      p <- draw_manhattanPlot(
        gwasFile = "../../data/results/gwasResult.json",
        gwasUrl = NULL,
        adj_method = "bonferroni",
        thresh_p = 0.05,
        filter_pAdj = 1,
        filter_nPoints = Inf,
        filter_quant = 1,
        chr = NA
      )
    })

    expect_true(all.equal(class(p), c("plotly", "htmlwidget")))

    expect_no_error({
      tmpF <- tempfile(fileext = ".html")
      p <- draw_manhattanPlot(
        gwasFile = "../../data/results/gwasResult.json",
        gwasUrl = NULL,
        adj_method = "bonferroni",
        thresh_p = 0.05,
        chr = NA,
        filter_pAdj = 1,
        filter_nPoints = Inf,
        filter_quant = 1,
        outFile = tmpF
      )
    })
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)

    expect_no_error({
      tmpF <- tempfile(fileext = ".html")
      p <- draw_manhattanPlot(
        gwasFile = "../../data/results/gwasResult.json",
        gwasUrl = NULL,
        adj_method = "bonferroni",
        thresh_p = 0.05,
        chr = "1",
        filter_pAdj = 1,
        filter_nPoints = Inf,
        filter_quant = 1,
        outFile = tmpF
      )
    })
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)

    expect_no_error({
      tmpF <- tempfile(fileext = ".png")
      p <- draw_manhattanPlot(
        gwasFile = "../../data/results/gwasResult.json",
        gwasUrl = NULL,
        adj_method = "bonferroni",
        thresh_p = 0.05,
        chr = NA,
        filter_pAdj = 1,
        filter_nPoints = Inf,
        filter_quant = 1,
        interactive = FALSE,
        outFile = tmpF
      )
    })
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)
  })



  # Test draw_manhattanPlot() with wrong parameters:
  goodParams <- list(
    gwasFile = "../../data/results/gwasResult.json",
    gwasUrl = NULL,
    adj_method = "bonferroni",
    thresh_p = 0.05,
    chr = NA,
    filter_pAdj = 1,
    filter_nPoints = Inf,
    filter_quant = 1,
    interactive = TRUE,
    outFile = NULL
  )

  wrongParamsL <- list(
    gwasFile = c("doNotExist"),
    gwasUrl = c("doNotExist"),
    adj_method = c("doNotExist"),
    thresh_p = list(-1, 1.02),
    chr = c("doNotExist", 50),
    filter_pAdj = c(-1, 2),
    filter_nPoints = -4,
    filter_quant = c(-1, 1.1),
    interactive = c("f"),
    outFile = list(c("f1", "f2"))
  )
  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i + 1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("draw_manhattanPlot, WrongParams", p, i, sep = "-")
      test_that(testName, {
        err <- expect_engineError({
          draw_manhattanPlot(
            gwasFile = params[["gwasFile"]],
            gwasUrl = params[["gwasUrl"]],
            adj_method = params[["adj_method"]],
            thresh_p = params[["thresh_p"]],
            chr = params[["chr"]],
            filter_pAdj = params[["filter_pAdj"]],
            filter_quant = params[["filter_quant"]],
            filter_nPoints = params[["filter_nPoints"]],
            interactive = params[["interactive"]],
            outFile = params[["outFile"]]
          )
        })
      })
    }
  }



  # adjust p-values ----
  test_that("GWAS Results adjust p-values", {
    expect_no_error({
      gwasAdjResults <- run_resAdjustment(
        gwasFile = "../../data/results/gwasResult.json",
        gwasUrl = NULL,
        adj_method = "bonferroni",
        filter_pAdj = 1,
        filter_nPoints = Inf,
        filter_quant = 1,
        outFile = tempfile()
      )
    })
    expect_true(class(gwasAdjResults$gwasAdjusted) == "json")
    expect_no_error({
      gwas <- readGWAS(gwasAdjResults$file)
    })
    expect_true(class(gwas) == "list")
    expect_true(all.equal(names(gwas), c("gwas", "metadata")))
    expect_true(class(gwas$gwas) == "data.frame")
    origGWAS <- readGWAS("../../data/results/gwasResult.json")
    expect_equal(nrow(gwas$gwas), nrow(origGWAS$gwas))
    expect_equal(ncol(gwas$gwas), ncol(origGWAS$gwas) + 1)
    expect_true(class(gwas$metadata) == "list")
    expect_true(
      all.equal(names(gwas$metadata), c(
        "genoFP",
        "phenoFP",
        "trait",
        "test",
        "fixed",
        "response",
        "thresh_maf",
        "thresh_callrate",
        "date",
        "adj_method"
      ))
    )
  })

  # Test run_resAdjustment() with wrong parameters:
  goodParams <- list(
    gwasFile = "../../data/results/gwasResult.json",
    gwasUrl = NULL,
    adj_method = "bonferroni",
    filter_pAdj = 1,
    filter_nPoints = Inf,
    filter_quant = 1,
    outFile = NULL
  )

  wrongParamsL <- list(
    gwasFile = c("doNotExist", NA),
    gwasUrl = c("doNotExist", NA),
    adj_method = c("doNotExist", NA),
    filter_pAdj = c(-1, 2),
    filter_nPoints = -4,
    filter_quant = c(-1, 1.1),
    outFile = list(c("f1", "f2"))
  )


  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i + 1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("run_resAdjustment, WrongParams", p, i, sep = "-")
      test_that(testName, {
        err <- expect_engineError({
          run_resAdjustment(
            gwasFile = params[["gwasFile"]],
            gwasUrl = params[["gwasUrl"]],
            adj_method = params[["adj_method"]],
            filter_pAdj = params[["filter_pAdj"]],
            filter_quant = params[["filter_quant"]],
            filter_nPoints = params[["filter_nPoints"]],
            outFile = params[["outFile"]]
          )
        })
      })
    }
  }
})
