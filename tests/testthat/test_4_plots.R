# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Unit test for LD plots

capture_output({

  # LDplot ----
  test_that("LD plot", {
    gDta <- readGenoData("../../data/geno/testMarkerData01.vcf.gz")
    expect_error({
      imgFile <- LDplot(gDta, 1, 20, file = NULL)
    }, NA)
    expect_error({
      imgFile <- LDplot(gDta, 42, 72, file = tempfile())
    }, NA)
    expect_error({
      imgFile <- LDplot(gDta, 1024, 1074, file = "testOutput/plot.png")
    }, NA)
    file.remove(imgFile)

    # errors
    expect_error({
      imgFile <- LDplot(gDta, 1.5, 20, file = NULL)
    },'`from` should be an positive integer of length 1' )
    expect_error({
      imgFile <- LDplot(gDta, 1, 20.1, file = NULL)
    },'`to` should be an positive integer of length 1' )
    expect_error({
      imgFile <- LDplot(gDta, "2", 20, file = NULL)
    }, '`from` should be an positive integer of length 1')
    expect_error({
      imgFile <- LDplot(gDta, 2, "20", file = NULL)
    }, '`to` should be an positive integer of length 1')

    expect_error({
      imgFile <- LDplot(gDta, 20, 2, file = NULL)
    }, '"from" should be lower than "to"')
    expect_error({
      imgFile <- LDplot(gDta, 20, 200, file = NULL)
    }, 'range size \\("to" - "from"\\) should be lower or equal than 50')
  })

  # manPlot ----
  files <- list(c(g = "../../data/geno/testMarkerData01.vcf.gz",
                  p = "../../data/pheno/testPhenoData01.csv"))
  for (file in files) {
    dta <- readData(file["g"], file["p"])
    for (trait in c("Flowering.time.at.Arkansas", "Awn.presence")) {
      for (test in c("score", "wald", "lrt")) {
        if (length(na.omit(unique(dta$phenoData[,trait]))) > 2) {
          resp <- "quantitative"
        } else {
          if (test != "score") {
            next
          }
          resp <- "binary"
        }
        resGwas <- gwas(dta,
                        trait,
                        test,
                        fixed = 0,
                        response = resp,
                        thresh_maf = 0.05,
                        thresh_callrate = 0.95)
        for (chr in c(NA, sample(unique(resGwas$chr), 2))) {
          for (interactive in c(TRUE, FALSE)) {
            testName <- paste("manPlot", resp, test, paste0("chr:", chr),
                              paste0('inter:', interactive), sep = "-")
            test_that(testName, {
              expect_error({
                mP <- manPlot(gwas = resGwas,
                              adj_method = "bonferroni",
                              thresh_p = 0.05,
                              chr = chr,
                              interactive = interactive,
                              filter_pAdj = 1,
                              filter_nPoints = Inf,
                              filter_quant = 1,
                              title = "test manPlot")
              }, NA)
              if (interactive) {
                expect_true(all.equal(class(mP), c("plotly", "htmlwidget")))
              } else {
                expect_true(is.null(mP))
              }
            })
          }
        }
      }
    }
  }

  # Test manPlot warnings: ----
  expect_warning({
  emptyResGwas <- gwas(dta,
                       trait,
                       test,
                       fixed = 0,
                       response = resp,
                       thresh_maf = 0.5,
                       thresh_callrate = 0.95)
  })
  test_that('manPlot warning, empty data', {
    expect_warning({
      mP <- manPlot(gwas = emptyResGwas,
                    adj_method = "bonferroni",
                    thresh_p = 0.05,
                    chr = NA,
                    interactive = TRUE,
                    filter_pAdj = 1,
                    filter_nPoints = Inf,
                    filter_quant = 0,
                    title = "test manPlot")
    }, 'There is no points to display')
  })


  for (inter in c(TRUE, FALSE)) {
    testName <- paste("manPlot warnings, interactive:", inter)
    test_that(testName, {
      expect_warning({
        mP <- manPlot(gwas = resGwas,
                      adj_method = "bonferroni",
                      thresh_p = 0.05,
                      chr = NA,
                      interactive = inter,
                      filter_pAdj = 10^-9,
                      title = "test manPlot")
      }, 'removed all the points of the graph')
      expect_warning({
        mP <- manPlot(gwas = resGwas,
                      adj_method = "bonferroni",
                      thresh_p = 0.05,
                      chr = NA,
                      interactive = inter,
                      filter_quant = 0,
                      title = "test manPlot")
      }, 'removed all the points of the graph')
      expect_warning({
        mP <- manPlot(gwas = resGwas,
                      adj_method = "bonferroni",
                      thresh_p = 0.05,
                      chr = NA,
                      interactive = inter,
                      filter_nPoints = 0,
                      title = "test manPlot")
      }, 'removed all the points of the graph')
    })
  }


  # Test manPlot() with wrong parameters: ----
  resGwas <- gwas(dta,
                  trait = "Flowering.time.at.Arkansas",
                  test = "score",
                  fixed = 0,
                  response = "quantitative",
                  thresh_maf = 0.05,
                  thresh_callrate = 0.95)
  goodParams <- list(gwas = resGwas,
                     adj_method = "bonferroni",
                     thresh_p = 0.05,
                     interactive = TRUE,
                     filter_pAdj = 1,
                     filter_nPoints = Inf,
                     filter_quant = 1,
                     chr = unique(resGwas$chr)[1])

  wrongParamsL <- list(p = list(c(runif(100, 0, 1), -0.1, 0.00015)),
                       adj_method = "do not exists",
                       thresh_p = list(-1, 1.02),
                       interactive = c('t'),
                       filter_pAdj = c(-1, 2),
                       filter_nPoints = -4,
                       filter_quant = c(-1, 1.1),
                       chr = c("doNotExists", 29))

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i + 1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("manPlot WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          manPlot(gwas = params[["gwas"]],
                  adj_method = params[["adj_method"]],
                  thresh_p = params[["thresh_p"]],
                  chr = params[["chr"]],
                  interactive = params[["interactive"]],
                  filter_pAdj = params[['filter_pAdj']],
                  filter_quant = params[['filter_quant']],
                  filter_nPoints = params[['filter_nPoints']])
        })
      })
    }
  }




  # pedNetwork ----
  pedFiles <- c('../../data/pedigree/testPedData_char.csv',
                '../../data/pedigree/testPedData_num.csv',
                '../../data/pedigree/testPedData_missFounder.csv')
  for (file in pedFiles) {
    test_that(paste("pedNetwork", basename(file)), {
      suppressWarnings({
        ped <- readPedData(file)
      })
      expect_error({
        pn <- pedNetwork(ped)
      }, NA)
      expect_true(all.equal(class(pn), c("forceNetwork", "htmlwidget")))
    })
  }

  for (file in pedFiles) {
    for (interactive in c(TRUE, FALSE)) {
      test_that(paste("relMatHeatmap - inter:", interactive, basename(file)), {
        suppressWarnings({
          ped <- readPedData(file)
          relMat <- pedRelMat(ped)
        })
        expect_error({
          hm <- relMatHeatmap(relMat, interactive)
        }, NA)
        if (interactive) {
          expect_true(all.equal(class(hm), c("plotly", "htmlwidget")))
        } else {
          expect_true(is.null(hm))
        }
      })
    }
  }


  # plotBlup ----
  progEstimFiles <- '../../data/results/progenyBlupEstimation.json'
  for (file in progEstimFiles) {
    test_that(paste("plotBlup", basename(file)), {
      progBlup <- readProgBlupEstim(file)
      expect_error({
        p <- plotBlup(progBlup)
      }, NA)
      expect_true(all.equal(class(p), c("plotly", "htmlwidget")))
      expect_type(p$x$visdat[[1]]()$Cross, 'character')

      for (s in c('asc', 'dec', '???')) {
        expect_error({
          p <- plotBlup(progBlup, sorting = s)
        }, NA)
        if (s != '???') {
          expect_s3_class(p$x$visdat[[1]]()$Cross, 'factor')
        } else {
          expect_type(p$x$visdat[[1]]()$Cross, 'character')
        }
      }

    })
  }



})
