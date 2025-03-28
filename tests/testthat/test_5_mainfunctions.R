# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test of main functions

capture.output({


  # GS model ----
  tests_cases <- list(
    test_01 = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                   pheno = "../../data/genomic_selection/pheno_train.csv",
                   trait = "pheno"),
    missing_genotype = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                            pheno = "../data/GS_pheno_train_missing_genotype.csv",
                            trait = "pheno"),
    missing_phenotype = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                             pheno = "../data/GS_pheno_train_missing_phenotype.csv",
                             trait = "pheno"),
    duplicated_snp_ids = list(geno = "../data/geno_G1_duplicated_ids.vcf.gz",
                              pheno = "../../data/genomic_selection/pheno_train.csv",
                              trait = "pheno"),
    unused_trait_badly_formatted = list(geno = "../data/geno_G1_duplicated_ids.vcf.gz",
                                        pheno = "../../data/genomic_selection/pheno_train_with_pheno2_trait_badly_formated.csv",
                                        trait = "pheno")
  )

  for (test in names(tests_cases)) {

    for (with_dominance in c(T,F)) {
      test_that(paste0("train_gs_model_main ", test, " dominance: ", with_dominance), {
        expect_no_error({
          model <- train_gs_model_main(genoFile = tests_cases[[test]]$geno,
                                       phenoFile = tests_cases[[test]]$pheno,
                                       genoUrl = NULL,
                                       phenoUrl = NULL,
                                       trait = tests_cases[[test]]$trait,
                                       with_dominance = with_dominance,
                                       thresh_maf = 0,
                                       outFile = tempfile(fileext = ".json"))
        })
        expect_equal(1,1) # testthat v3.2.2 have a bug with `expect_no_error`
                          # if it succeed the test is considered as "skipped".

      })
    }
  }


  tests_bad_cases <- list(
    no_heterozygot = list(geno = "../../data/geno/testMarkerData01.vcf.gz",
                          pheno = "../../data/pheno/testPhenoData01.csv",
                          trait = "Flowering.time.at.Arkansas")
  )
  for (test in names(tests_bad_cases)) {

    test_that(paste0("train_gs_model_main Bad case ", test), {
      expect_engineError({
        model <- train_gs_model_main(genoFile = tests_bad_cases[[test]]$geno,
                                     phenoFile = tests_bad_cases[[test]]$pheno,
                                     genoUrl = NULL,
                                     phenoUrl = NULL,
                                     trait = tests_bad_cases[[test]]$trait,
                                     with_dominance = TRUE,
                                     thresh_maf = 0,
                                     outFile = tempfile(fileext = ".json"))
      })
    })
  }


  # GS cross-validation ----
  tests_cases <- list(
    test_01 = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                   pheno = "../../data/genomic_selection/pheno_train.csv",
                   trait = "pheno"),
    missing_genotype = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                            pheno = "../data/GS_pheno_train_missing_genotype.csv",
                            trait = "pheno"),
    missing_phenotype = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                             pheno = "../data/GS_pheno_train_missing_phenotype.csv",
                             trait = "pheno")
  )

  for (test in names(tests_cases)) {

    test_that(paste0("cross_validation_evaluation_main ", test), {
      for (with_dominance in c(T,F)) {
        expect_no_error({
          evaluation <- cross_validation_evaluation_main(genoFile = tests_cases[[test]]$geno,
                                                         phenoFile = tests_cases[[test]]$pheno,
                                                         genoUrl = NULL,
                                                         phenoUrl = NULL,
                                                         trait = tests_cases[[test]]$trait,
                                                         with_dominance = with_dominance,
                                                         n_folds = 10,
                                                         n_repetitions = 5,
                                                         thresh_maf = 0,
                                                         outFile = tempfile(fileext = ".json")
          )
        })
      }
    })
  }


  # GS evaluation plot ----
  tests_cases <- list(
    test_01 = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                   pheno = "../../data/genomic_selection/pheno_train.csv",
                   trait = "pheno"),
    missing_genotype = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                            pheno = "../data/GS_pheno_train_missing_genotype.csv",
                            trait = "pheno"),
    missing_phenotype = list(geno = "../../data/genomic_selection/geno_G1.vcf.gz",
                             pheno = "../data/GS_pheno_train_missing_phenotype.csv",
                             trait = "pheno")
  )

  for (test in names(tests_cases)) {
    test_that(paste0("cross_validation_evaluation_main ", test), {
      for (with_dominance in c(T,F)) {
        evaluation <- cross_validation_evaluation_main(genoFile = tests_cases[[test]]$geno,
                                                       phenoFile = tests_cases[[test]]$pheno,
                                                       trait = tests_cases[[test]]$trait,
                                                       with_dominance = with_dominance,
                                                       n_folds = 10,
                                                       n_repetitions = 5,
                                                       thresh_maf = 0,
                                                       outFile = tempfile(fileext = ".json"))
        expect_no_error({
          draw_evaluation_plot(evaluation$file,
                               tempfile(fileext = ".html"))
        })
      }

    })
  }

  # GS predictions ----

  tests_cases <- list(
    additive = list(geno = "../../data/genomic_selection/geno_G2.vcf.gz",
                   markerEffects = "../../data/results/GS_model_additive.json"),
    dominance = list(geno = "../../data/genomic_selection/geno_G2.vcf.gz",
                     markerEffects = "../../data/results/GS_model_dominance.json"),
    several_traits = list(geno = "../../data/geno/breedGame_phasedGeno.vcf.gz",
                          markerEffects = "../../data/markerEffects/breedGame_markerEffects_2traits.json")
  )

  for (test in names(tests_cases)) {
    test_that(paste0("train_gs_model_main ", test), {
      expect_no_error({
        predictions <- predict_gs_model_main(genoFile = tests_cases[[test]]$geno,
                                             genoUrl = NULL,
                                             markerEffectsFile = tests_cases[[test]]$markerEffects,
                                             markerEffectsUrl = NULL,
                                             outFile = tempfile(fileext = ".csv")
        )
      })
    })
  }





  # run GWAS ----
  files <- list(
    test_01 = list(geno =  "../../data/geno/testMarkerData01.vcf.gz" ,
                   pheno = "../../data/pheno/testPhenoData01.csv",
                   trait = "Flowering.time.at.Arkansas"),
    test_singleTrait = list(geno = "../data/Result_genos_hd_subset-initialColl.vcf.gz",
                            pheno = "../data/resistance_initColl.csv",
                            trait = "resist"),
    missingSNPid = list(geno = "../data/Result_genos_hd_subset-initialColl_missingId.vcf.gz",
                        pheno = "../data/resistance_initColl.csv",
                        trait = "resist"),
    unused_trait_badly_formatted = list(geno = "../data/geno_G1_duplicated_ids.vcf.gz",
                                        pheno = "../../data/genomic_selection/pheno_train_with_pheno2_trait_badly_formated.csv",
                                        trait = "pheno")
  )
  for (test in names(files)) {

    test_that(paste0("Run GWAS", test), {
      expect_no_error({
        gwas_results <- run_gwas(genoFile = files[[test]]$geno,
                                 phenoFile = files[[test]]$pheno,
                                 genoUrl = NULL,
                                 phenoUrl = NULL,
                                 trait = files[[test]]$trait,
                                 test = "score",
                                 fixed = 0,
                                 response = "quantitative",
                                 thresh_maf = 0.05,
                                 thresh_callrate = 0.95,
                                 outFile = tempfile(fileext = ".json"))
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
  }


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
                       phenoFile = c("/pheno-do-not-exist", NA,
                                     "../data/pheno_duplicated.csv",
                                     "../data/pheno_no_common_ind_with_geno.csv",
                                     "../data/pheno_inconsistent_both_numeric_and_character_values.csv"),
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
        err <- expect_engineError({
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


  # draw_manhattanPlot ----
  test_that("Draw Manhattan Plot", {
    expect_no_error({
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
                              chr = NA)
    })

    expect_true(all.equal(class(p), c("plotly", "htmlwidget")))

    expect_no_error({
      tmpF <- tempfile(fileext = ".html")
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = NA,
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
                              outFile = tmpF)
    })
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)

    expect_no_error({
      tmpF <- tempfile(fileext = ".html")
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = "1",
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
                              outFile = tmpF)
    })
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)

    expect_no_error({
      tmpF <- tempfile(fileext = ".png")
      p <- draw_manhattanPlot(gwasFile = "../../data/results/gwasResult.json",
                              gwasUrl = NULL,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = NA,
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
                              interactive = FALSE,
                              outFile = tmpF)
    })
    expect_true(file.exists(tmpF))
    expect_true(file.info(tmpF)$size > 0)
  })



  # Test draw_manhattanPlot() with wrong parameters:
  goodParams <- list(gwasFile = "../../data/results/gwasResult.json",
                     gwasUrl = NULL,
                     adj_method = "bonferroni",
                     thresh_p = 0.05,
                     chr = NA,
                     filter_pAdj = 1,
                     filter_nPoints = Inf,
                     filter_quant = 1,
                     interactive = TRUE,
                     outFile = NULL)

  wrongParamsL <- list(gwasFile = c("doNotExist"),
                       gwasUrl = c("doNotExist"),
                       adj_method = c("doNotExist"),
                       thresh_p = list(-1, 1.02),
                       chr = c("doNotExist", 50),
                       filter_pAdj = c(-1, 2),
                       filter_nPoints = -4,
                       filter_quant = c(-1, 1.1),
                       interactive = c('f'),
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
        err <- expect_engineError({
          draw_manhattanPlot(gwasFile = params[["gwasFile"]],
                             gwasUrl = params[["gwasUrl"]],
                             adj_method = params[["adj_method"]],
                             thresh_p = params[["thresh_p"]],
                             chr = params[["chr"]],
                             filter_pAdj = params[['filter_pAdj']],
                             filter_quant = params[['filter_quant']],
                             filter_nPoints = params[['filter_nPoints']],
                             interactive = params[["interactive"]],
                             outFile = params[["outFile"]])
        })
      })
    }
  }

  # Draw LD plot ----
  test_that("Draw LD plot", {
    expect_no_error({
      imgFile <- draw_ldPlot(genoFile = "../../data/geno/testMarkerData01.vcf.gz",
                             genoUrl = NULL,
                             from = 42,
                             to = 62,
                             outFile = tempfile(fileext = ".png"))
    })
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
        err <- expect_engineError({
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
    expect_no_error({
      gwasAdjResults <- run_resAdjustment(gwasFile = "../../data/results/gwasResult.json",
                                          gwasUrl = NULL,
                                          adj_method = "bonferroni",
                                          filter_pAdj = 1,
                                          filter_nPoints = Inf,
                                          filter_quant = 1,
                                          outFile = tempfile())
    })
    expect_true(class(gwasAdjResults$gwasAdjusted) == "json")
    expect_no_error({
      gwas <- readGWAS(gwasAdjResults$file)
    })
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
                     filter_pAdj = 1,
                     filter_nPoints = Inf,
                     filter_quant = 1,
                     outFile = NULL)

  wrongParamsL <- list(gwasFile = c("doNotExist", NA),
                       gwasUrl = c("doNotExist", NA),
                       adj_method = c("doNotExist", NA),
                       filter_pAdj = c(-1, 2),
                       filter_nPoints = -4,
                       filter_quant = c(-1, 1.1),
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
        err <- expect_engineError({
          run_resAdjustment(gwasFile = params[["gwasFile"]],
                            gwasUrl = params[["gwasUrl"]],
                            adj_method = params[["adj_method"]],
                            filter_pAdj = params[['filter_pAdj']],
                            filter_quant = params[['filter_quant']],
                            filter_nPoints = params[['filter_nPoints']],
                            outFile = params[["outFile"]])
        })
      })
    }
  }





  # calc_pedRelMat ----
  pedFiles <- c('../../data/pedigree/testPedData_char.csv',
                '../../data/pedigree/testPedData_num.csv',
                '../../data/pedigree/testPedData_missFounder.csv')
  formats <- c('csv', 'json')
  for (file in pedFiles) {
    for (format in formats) {
      test_that(paste("calc_pedRelMat", format, basename(file)), {
        expect_no_error({
          suppressWarnings({
            relMat_results <- calc_pedRelMat(pedFile = file,
                                             pedUrl = NULL,
                                             header = TRUE,
                                             unknown_string = '',
                                             outFormat = format)
          })
        })

        expect_true(class(relMat_results) == "list")
        expect_identical(names(relMat_results),
                         c("relMat", "metadata", "file"))
        expect_true(is.matrix(relMat_results$relMat))
        expect_identical(names(relMat_results$metadata),
                         c("info", "date", "nInds", "pedFP"))
        expect_true(file.exists(relMat_results$file))
      })
    }
  }

  # calc_genoRelMat ----
  genoFiles <- c('../../data/geno/breedGame_geno.vcf.gz')
  formats <- c('csv', 'json')
  for (format in formats) {
    for (file in genoFiles) {
      test_that(paste("calc_genoRelMat", format, basename(file)), {
        expect_no_error({
          relMat_results <- calc_genoRelMat(genoFile = file,
                                            genoUrl = NULL,
                                            outFormat = format)
        })

        expect_true(class(relMat_results) == "list")
        expect_identical(names(relMat_results),
                         c("relMat", "metadata", "file"))
        expect_true(is.matrix(relMat_results$relMat))
        expect_identical(names(relMat_results$metadata),
                         c("info", "date", "nInds", "genoFP"))
        expect_true(file.exists(relMat_results$file))
      })
    }
  }



  # calc_combinedRelMat ----
  rm_files <- list(
    breedGameFiles = list(
      ped_rm = c('../../data/results/breedGame_pedRelMat.csv'),
      geno_rm = c('../../data/results/breedGame_genoRelMat.csv')
    )
  )
  formats <- c('csv', 'json')
  methods <- c('Legarra', 'Martini')
  tau <- list(Legarra = NULL, Martini = 0.4)
  omega <- list(Legarra = NULL, Martini = 0.7)
  for (format in formats) {
    for (method in methods) {
      for (filesName in names(rm_files)) {
        test_that(paste("calc_genoRelMat", format, method, filesName), {
          files <- rm_files[[filesName]]
          expect_no_error({
            relMat_results <- calc_combinedRelMat(
              pedRelMatFile = files$ped_rm,
              genoRelMatFile = files$geno_rm,
              method = method,
              tau = tau[[method]],
              omega = omega[[method]],
              outFormat = format
            )
          })

          expect_true(class(relMat_results) == "list")
          expect_identical(names(relMat_results),
                           c("relMat", "metadata", "file"))
          expect_true(is.matrix(relMat_results$relMat))
          expect_identical(names(relMat_results$metadata),
                           c("info", "date", "nInds",
                             "geno_relMatFP" , "ped_relMatFP"))
          expect_true(file.exists(relMat_results$file))
        })
      }
    }
  }





  # draw_pedNetrowk ----
  for (file in pedFiles) {
    test_that(paste("draw_pedNetwork", basename(file)), {
      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        suppressWarnings({
          pedNet <- draw_pedNetwork(pedFile = file,
                                    pedUrl = NULL,
                                    header = TRUE,
                                    unknown_string = '',
                                    outFile = tmpF)
        })
      })
      expect_identical(class(pedNet),
                       c("forceNetwork", "htmlwidget"))
      expect_true(file.exists(tmpF))
      expect_true(file.info(tmpF)$size > 0)
    })
  }

  # draw_relMat ----
  relMatFiles <- c('../../data/results/pedigreeRelationship.csv',
                   '../../data/results/pedigreeRelationship.json')
  for (file in relMatFiles) {
    for (inter in c(TRUE, FALSE)) {

      test_that(paste("draw_relHeatmap inter:", inter, basename(file)), {
        ext <- ".html"
        if (!inter) {
          ext <- ".png"
        }
        tmpF <- tempfile(fileext = ext)
        expect_no_error({
          pedNet <- draw_relHeatmap(relMatFile = file,
                                    relMatUrl = NULL,
                                    interactive = inter,
                                    outFile = tmpF)
        })
        expect_true(file.exists(tmpF))
        expect_true(file.info(tmpF)$size > 0)
      })
    }
  }



  # crossingSimulation ----
  crossigSimDataList <- list(
    completeDta = list(
      phasedGenoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      SNPcoordFile = '../../data/SNPcoordinates/breedingGame_SNPcoord.csv',
      crossTableFile = '../../data/crossingTable/breedGame_crossTable.csv',
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingChr = list(
      phasedGenoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingChr.csv',
      crossTableFile = '../../data/crossingTable/breedGame_crossTable.csv',
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingSNPid = list(
      phasedGenoFile = '../data/breedGame_phasedGeno_missingIds.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingIds.csv',
      crossTableFile = '../../data/crossingTable/breedGame_crossTable.csv',
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingPhysPos = list(
      phasedGenoFile = '../data/breedGame_phasedGeno_missingPos.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingPhysPos.csv',
      crossTableFile = '../../data/crossingTable/breedGame_crossTable.csv',
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    ),
    missingPhysPosColumns = list(
      phasedGenoFile = '../data/breedGame_phasedGeno_missingPos.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingPhysPos_columns.csv',
      crossTableFile = '../../data/crossingTable/breedGame_crossTable.csv',
      nCross = 3,
      outFile = tempfile(fileext = ".vcf.gz")
    )
  )

  for (i in seq_along(crossigSimDataList)) {
    data <- crossigSimDataList[[i]]
    phasedGenoFile <- data$phasedGenoFile
    SNPcoordFile <- data$SNPcoordFile
    crossTableFile <- data$crossTableFile
    nCross <- data$nCross
    outFile <- data$outFile

    test_that(paste('crossingSimulation -', names(crossigSimDataList)[i]), {
      if (grepl(pattern = "missingPhysPos|missingChr", names(crossigSimDataList)[i])) {
        expect_warning({
          createdFile <- crossingSimulation(genoFile = phasedGenoFile,
                                            crossTableFile = crossTableFile,
                                            SNPcoordFile = SNPcoordFile,
                                            nCross = nCross,
                                            outFile = outFile)
        })
      } else {
        expect_no_error({
          createdFile <- crossingSimulation(genoFile = phasedGenoFile,
                                            crossTableFile = crossTableFile,
                                            SNPcoordFile = SNPcoordFile,
                                            nCross = nCross,
                                            outFile = outFile)
        })
      }
      expect_equal(createdFile, outFile)
      expect_no_error({
        createdGeno <- readGenoData(createdFile)
      })
      expect_no_error({
        createdPhasedGeno <- readPhasedGeno(createdFile)
      })

      # individual names
      crossTable <- readCrossTable(crossTableFile)
      expect_equal(nrow(createdGeno@ped), nrow(crossTable) * nCross)
      createdIndBaseName <- unique(gsub('-[0-9]*$', '', createdGeno@ped$id))
      expect_equal(sort(createdIndBaseName), sort(crossTable$names))

      # SNPs: check original and new vcf files' snp position are the same
      g <- readPhasedGeno(phasedGenoFile)
      snpPhysPos_newVcf <- createdPhasedGeno$SNPcoord[order(createdPhasedGeno$SNPcoord$SNPid),]
      snpPhysPos_ref <- g$SNPcoord[order(g$SNPcoord$SNPid),]
      expect_equal(snpPhysPos_newVcf$SNPid, snpPhysPos_ref$SNPid)
      expect_equal(snpPhysPos_newVcf$physPos, snpPhysPos_ref$physPos)

      if (grepl('missingPhysPos', names(crossigSimDataList)[i])) {
        expect_true(all(is.na(createdPhasedGeno$SNPcoord$physPos)))
      }
    })
  }

  test_that('crossingSimulation with downloaded files', {
    as_url <- function(file) {
      file <- normalizePath(file)
      file <- paste0("file://", file)
      file
    }
    expect_warning({
      createdFileWithoutChrInfo <- crossingSimulation(
        genoUrl = as_url(phasedGenoFile),
        crossTableUrl = as_url(crossTableFile),
        SNPcoordUrl = as_url(SNPcoordFile),
        nCross = nCross,
        outFile = outFile)
    }, 'output "file" already exists. This file will be overwritten.' )
  })
  unlink(outFile)

  inconsistentSNPFile <- '../data/inconsistent_SNPcoord_2.csv'
  test_that('crossingSimulation inconsistent SNPs', {
    err <- expect_engineError({
      createdFile <- crossingSimulation(
        genoFile = phasedGenoFile,
        crossTableFile = crossTableFile,
        SNPcoordFile = inconsistentSNPFile,
        nCross = nCross,
        outFile = outFile)
    })
  })



  # calc_progenyBlupEstimation ----
  progBlupDataList <- list(
    completeDta = list(
      genoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      SNPcoordFile = '../../data/SNPcoordinates/breedingGame_SNPcoord.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects.csv',
      outFile = tempfile(fileext = ".json")
    ),
    missingChr = list(
      genoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingChr.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects.csv',
      outFile = tempfile(fileext = ".json")
    ),
    missingSNPid = list(
      genoFile = '../data/breedGame_phasedGeno_missingIds.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingIds.csv',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      markerEffectsFile = '../data/breedGame_markerEffects_imputedSNPid.csv',
      outFile = tempfile(fileext = ".json")
    ),
    missingPhysPos = list(
      genoFile = '../data/breedGame_phasedGeno_missingPos.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingPhysPos.csv',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects.csv',
      outFile = tempfile(fileext = ".json")
    ),
    missingPhysPosColumns = list(
      genoFile = '../data/breedGame_phasedGeno_missingPos.vcf.gz',
      SNPcoordFile = '../data/breedingGame_SNPcoord_missingPhysPos_columns.csv',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects.csv',
      outFile = tempfile(fileext = ".json")
    ),
    completeDta_severalTraits = list(
      genoFile = '../../data/geno/breedGame_phasedGeno.vcf.gz',
      crossTableFile = '../../data/crossingTable/breedGame_small_crossTable.csv',
      SNPcoordFile = '../../data/SNPcoordinates/breedingGame_SNPcoord.csv',
      markerEffectsFile = '../../data/markerEffects/breedGame_markerEffects_2traits.json',
      outFile = tempfile(fileext = ".json")
    )
  )

  for (i in seq_along(progBlupDataList)) {

    genoFile <- progBlupDataList[[i]]$genoFile
    crossTableFile <-  progBlupDataList[[i]]$crossTableFile
    SNPcoordFile <-  progBlupDataList[[i]]$SNPcoordFile
    markerEffectsFile <-  progBlupDataList[[i]]$markerEffectsFile
    outFile <-  progBlupDataList[[i]]$outFile

    test_that(paste('calc_progenyBlupEstimation -', names(progBlupDataList)[i]), {
      if (grepl(pattern = "missingPhysPos|missingChr", names(crossigSimDataList)[i])) {
        expect_warning({
          projBlups_list <- calc_progenyBlupEstimation(
            genoFile = genoFile,
            crossTableFile = crossTableFile,
            SNPcoordFile = SNPcoordFile,
            markerEffectsFile = markerEffectsFile,
            outFile = outFile
          )
        })
      } else {
        expect_no_error({
          projBlups_list <- calc_progenyBlupEstimation(
            genoFile = genoFile,
            crossTableFile = crossTableFile,
            SNPcoordFile = SNPcoordFile,
            markerEffectsFile = markerEffectsFile,
            outFile = outFile
          )
        })
      }
      expect_is(projBlups_list, "list")
      lapply(projBlups_list, function(projBlups){
        expect_is(projBlups, "list")
        expect_identical(names(projBlups),
                         c("ind1", "ind2", "blup_exp", "blup_var", "cov"))
      })
    })

  }


  # draw_progBlupsPlot ----
  progEstimFiles <- c('../../data/results/progenyBlupEstimation.json',
                      '../data/progenyBlupEstimation_oldVersion.json')
  for (file in progEstimFiles) {
    test_that(paste("draw_progBlupsPlot", basename(file)), {

      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        p1 <- draw_progBlupsPlot(progEstimFile = file,
                                 outFile = tmpF)
        unlink(tmpF)
      })

      for (s in c('asc', 'dec')) {
        expect_no_error({
          p1 <- draw_progBlupsPlot(progEstimFile = file,
                                   sorting = s,
                                   outFile = tmpF)
          unlink(tmpF)
        })
      }
      expect_warning({
        p1 <- draw_progBlupsPlot(progEstimFile = file,
                                 sorting = "??????",
                                 outFile = tmpF)
        unlink(tmpF)
      }, "Sorting method not recognised, x axis will be sorted alphabetically. Recognised methods are: alpha, asc, dec")
    })
  }

  progEstimFiles_severalTraits <- '../../data/results/progenyBlupEstimation_2traits.json'
  for (file in progEstimFiles_severalTraits) {
    test_that(paste("draw_progBlupsPlot several traits", basename(file)), {

      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        p1 <- draw_progBlupsPlot(progEstimFile = file,
                                 trait = 'trait1',
                                 outFile = tmpF)
        unlink(tmpF)
      })

      for (s in c('asc', 'dec')) {
        expect_no_error({
          p1 <- draw_progBlupsPlot(progEstimFile = file,
                                   sorting = s,
                                   trait = 'trait1',
                                   outFile = tmpF)
          unlink(tmpF)
        })
      }
      expect_warning({
        p1 <- draw_progBlupsPlot(progEstimFile = file,
                                 sorting = "??????",
                                 trait = 'trait1',
                                 outFile = tmpF)
        unlink(tmpF)
      }, "Sorting method not recognised, x axis will be sorted alphabetically. Recognised methods are: alpha, asc, dec")
    })
  }


  # draw_progBlupsPlot_2traits ----
  for (file in progEstimFiles_severalTraits) {
    test_that(paste("draw_progBlupsPlot_2traits several traits", basename(file)), {

      tmpF <- tempfile(fileext = ".html")
      expect_no_error({
        p1 <- draw_progBlupsPlot_2traits(progEstimFile = file,
                                         x_trait = 'trait1',
                                         y_trait = 'trait2',
                                         confidenceLevel = 0.95,
                                         x_suffix = "",
                                         y_suffix = "",
                                         ellipses_npoints = 100,
                                         outFile = tmpF)
        unlink(tmpF)
      })
    })
  }

})
