# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test of main functions

capture.output({

  # run GWAS ----
  files <- list(
    test_01 = list(geno =  "../../data/geno/testMarkerData01.vcf.gz" ,
                   pheno = "../../data/pheno/testPhenoData01.csv"),
    test_singleTrait = list(geno = "../data/Result_genos_hd_subset-initialColl.vcf.gz",
                          pheno = "../data/resistance_initColl.csv")
  )
  for (test in names(files)) {

    test_that(paste0("Run GWAS", test), {
      if (test == 'test_01') {
        trait <- "Flowering.time.at.Arkansas"
      } else if (test == 'test_singleTrait') {
        trait <- "resist"
      }
      expect_error({
        gwas_results <- run_gwas(genoFile = files[[test]]$geno,
                                 phenoFile = files[[test]]$pheno,
                                 genoUrl = NULL,
                                 phenoUrl = NULL,
                                 trait = trait,
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
                       phenoFile = c("/pheno-do-not-exist", NA, "../data/pheno_duplicated.csv"),
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


  # draw_manhattanPlot ----
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
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
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
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
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
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
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
                     filter_pAdj = 1,
                     filter_nPoints = Inf,
                     filter_quant = 1,
                     interactive = TRUE,
                     outFile = NULL)

  wrongParamsL <- list(gwasFile = c("doNotExist", NA),
                       gwasUrl = c("doNotExist", NA),
                       adj_method = c("doNotExist", NA),
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
        expect_error({
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
                                          filter_pAdj = 1,
                                          filter_nPoints = Inf,
                                          filter_quant = 1,
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
        expect_error({
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
      test_that(paste("calc_pedRelMAt", format, basename(file)), {
        expect_error({
          suppressWarnings({
          relMat_results <- calc_pedRelMAt(pedFile = file,
                                           pedUrl = NULL,
                                           header = TRUE,
                                           unknown_string = '',
                                           outFormat = format)
          })
        }, NA)

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
      test_that(paste("calc_genoRelMAt", format, basename(file)), {
        expect_error({
          relMat_results <- calc_genoRelMAt(genoFile = file,
                                            genoUrl = NULL,
                                            outFormat = format)
        }, NA)

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
        test_that(paste("calc_genoRelMAt", format, method, filesName), {
          files <- rm_files[[filesName]]
          expect_error({
            relMat_results <- calc_combinedRelMat(
              pedRelMatFile = files$ped_rm,
              genoRelMatFile = files$geno_rm,
              method = method,
              tau = tau[[method]],
              omega = omega[[method]],
              outFormat = format
            )
          }, NA)

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
      expect_error({
        suppressWarnings({
        pedNet <- draw_pedNetwork(pedFile = file,
                                  pedUrl = NULL,
                                  header = TRUE,
                                  unknown_string = '',
                                  outFile = tmpF)
        })
      }, NA)
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
        tmpF <- tempfile()
        expect_error({
          pedNet <- draw_relHeatmap(relMatFile = file,
                                    relMatUrl = NULL,
                                    interactive = inter,
                                    outFile = tmpF)
        }, NA)
        expect_true(file.exists(tmpF))
        expect_true(file.info(tmpF)$size > 0)
      })
    }
  }



  # crossingSimulation ----
    phasedGenoFile <- '../../data/geno/breedGame_phasedGeno.vcf.gz'
    SNPcoordFile <- '../../data/SNPcoordinates/breedingGame_SNPcoord.csv'
    crossTableFile <- '../../data/crossingTable/breedGame_crossTable.csv'
    chrInfoFile <- '../../data/chromosomesInformation/breedingGame_chrInfo.csv'
    nCross <- 3
    outFile <- tempfile(fileext = ".vcf.gz")

  test_that('crossingSimulation', {
    expect_error({
      createdFile <- crossingSimulation(genoFile = phasedGenoFile,
                                        crossTableFile = crossTableFile,
                                        SNPcoordFile = SNPcoordFile,
                                        chrInfoFile = chrInfoFile,
                                        nCross = nCross,
                                        outFile = outFile)
    }, NA)
    expect_equal(createdFile, outFile)
    expect_error({
      createdGeno <- readGenoData(createdFile)
    }, NA)
    expect_error({
      createdPhasedGeno <- readPhasedGeno(createdFile)
    }, NA)

    # individual names
    crossTable <- readCrossTable(crossTableFile)
    expect_equal(nrow(createdGeno@ped), nrow(crossTable) * nCross)
    createdIndBaseName <- unique(gsub('-[0-9]*$', '', createdGeno@ped$id))
    expect_equal(sort(createdIndBaseName), sort(crossTable$names))

    # SNPs
    SNPcoord <- readSNPcoord(SNPcoordFile)
    expect_equal(sort(createdPhasedGeno$SNPcoord$SNPid), sort(SNPcoord$SNPid))
    expect_equal(createdPhasedGeno$SNPcoord$physPos, SNPcoord$physPos)

  })

  test_that('crossingSimulation with downloaded files', {
    as_url <- function(file) {
      file <- normalizePath(file)
      file <- paste0("file://", file)
      file
    }
    expect_error({
      createdFileWithoutChrInfo <- crossingSimulation(
        genoUrl = as_url(phasedGenoFile),
        crossTableUrl = as_url(crossTableFile),
        SNPcoordUrl = as_url(SNPcoordFile),
        chrInfoUrl = as_url(chrInfoFile),
        nCross = nCross,
        outFile = outFile)
    }, NA)
  })

  test_that('crossingSimulation no chrInfo', {
    expect_error({
      createdFileWithoutChrInfo <- crossingSimulation(
        genoFile = phasedGenoFile,
        crossTableFile = crossTableFile,
        SNPcoordFile = SNPcoordFile,
        chrInfoFile = NULL,
        nCross = nCross,
        outFile = outFile)
    }, NA)
  })

  inconsistentSNPFile <- '../data/inconsistent_SNPcoord_2.csv'
  test_that('crossingSimulation inconsistent SNPs', {
    expect_error({
      createdFile <- crossingSimulation(
        genoFile = phasedGenoFile,
        crossTableFile = crossTableFile,
        SNPcoordFile = inconsistentSNPFile,
        chrInfoFile = chrInfoFile,
        nCross = nCross,
        outFile = outFile)
    }, paste("SNP's position order should be similar when sorted by",
             "physical position and by linkage map position."))
  })

})
