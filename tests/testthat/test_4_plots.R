# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Unit test for LD plots

capture_output({
  # LDplot ----
  test_that("LD plot", {
    gDta <- readGenoData("../../data/geno/testMarkerData01.vcf.gz")
    expect_no_error({
      imgFile <- LDplot(gDta, 1, 20, file = NULL)
    })
    expect_no_error({
      imgFile <- LDplot(gDta, 42, 72, file = tempfile())
    })
    expect_no_error({
      imgFile <- LDplot(gDta, 1024, 1074, file = tempfile(fileext = ".png"))
    })
    file.remove(imgFile)
  })

  # manPlot ----
  files <- list(c(
    g = "../../data/geno/testMarkerData01.vcf.gz",
    p = "../../data/pheno/testPhenoData01.csv"
  ))
  for (file in files) {
    dta <- readData(file["g"], file["p"])
    for (trait in c("Flowering.time.at.Arkansas", "Awn.presence")) {
      for (test in c("score", "wald", "lrt")) {
        if (length(na.omit(unique(dta$phenoData[, trait]))) > 2) {
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
          thresh_callrate = 0.95
        )
        for (chr in c(NA, sample(unique(resGwas$chr), 2))) {
          for (interactive in c(TRUE, FALSE)) {
            testName <- paste("manPlot", resp, test, paste0("chr:", chr),
              paste0("inter:", interactive),
              sep = "-"
            )
            test_that(testName, {
              expect_no_error({
                mP <- manPlot(
                  gwas = resGwas,
                  adj_method = "bonferroni",
                  thresh_p = 0.05,
                  chr = chr,
                  interactive = interactive,
                  filter_pAdj = 1,
                  filter_nPoints = Inf,
                  filter_quant = 1,
                  title = "test manPlot"
                )
              })
              if (interactive) {
                expect_true(all.equal(class(mP), c("plotly", "htmlwidget")))
              } else {
                expect_true(is.null(mP))
              }
            })

            test_that(paste(testName, "with some missing p values"), {
              resGwas$p[c(2, 4, 10)] <- NA
              expect_no_error({
                mP <- manPlot(
                  gwas = resGwas,
                  adj_method = "bonferroni",
                  thresh_p = 0.05,
                  chr = chr,
                  interactive = interactive,
                  filter_pAdj = 1,
                  filter_nPoints = Inf,
                  filter_quant = 1,
                  title = "test manPlot"
                )
              })
              if (interactive) {
                expect_true(all.equal(class(mP), c("plotly", "htmlwidget")))
              } else {
                expect_true(is.null(mP))
              }
            })

            test_that(paste(testName, "with all missing p values"), {
              resGwas$p <- NA
              expect_warning(
                {
                  mP <- manPlot(
                    gwas = resGwas,
                    adj_method = "bonferroni",
                    thresh_p = 0.05,
                    chr = chr,
                    interactive = interactive,
                    filter_pAdj = 1,
                    filter_nPoints = Inf,
                    filter_quant = 1,
                    title = "test manPlot"
                  )
                },
                regexp = "There is no points to display"
              )
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
      thresh_callrate = 0.95
    )
  })
  test_that("manPlot warning, empty data", {
    expect_warning(
      {
        mP <- manPlot(
          gwas = emptyResGwas,
          adj_method = "bonferroni",
          thresh_p = 0.05,
          chr = NA,
          interactive = TRUE,
          filter_pAdj = 1,
          filter_nPoints = Inf,
          filter_quant = 0,
          title = "test manPlot"
        )
      },
      "There is no points to display"
    )
  })


  for (inter in c(TRUE, FALSE)) {
    testName <- paste("manPlot warnings, interactive:", inter)
    test_that(testName, {
      expect_warning(
        {
          mP <- manPlot(
            gwas = resGwas,
            adj_method = "bonferroni",
            thresh_p = 0.05,
            chr = NA,
            interactive = inter,
            filter_pAdj = 10^-9,
            title = "test manPlot"
          )
        },
        "removed all the points of the graph"
      )
      expect_warning(
        {
          mP <- manPlot(
            gwas = resGwas,
            adj_method = "bonferroni",
            thresh_p = 0.05,
            chr = NA,
            interactive = inter,
            filter_quant = 0,
            title = "test manPlot"
          )
        },
        "removed all the points of the graph"
      )
      expect_warning(
        {
          mP <- manPlot(
            gwas = resGwas,
            adj_method = "bonferroni",
            thresh_p = 0.05,
            chr = NA,
            interactive = inter,
            filter_nPoints = 0,
            title = "test manPlot"
          )
        },
        "removed all the points of the graph"
      )
    })
  }

  # pedNetwork ----
  pedFiles <- c(
    "../../data/pedigree/testPedData_char.csv",
    "../../data/pedigree/testPedData_num.csv",
    "../../data/pedigree/testPedData_missFounder.csv"
  )
  for (file in pedFiles) {
    test_that(paste("pedNetwork", basename(file)), {
      suppressWarnings({
        ped <- readPedData(file)
      })
      expect_no_error({
        pn <- pedNetwork(ped)
      })
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
        expect_no_error({
          hm <- relMatHeatmap(relMat, interactive)
        })
        if (interactive) {
          expect_true(all.equal(class(hm), c("plotly", "htmlwidget")))
        } else {
          expect_true(is.null(hm))
        }
      })
    }
  }


  # plotBlup_1trait ----
  progEstimFiles <- c(
    "../../data/results/progenyBlupEstimation.json",
    "../data/progenyBlupEstimation_oldVersion.json"
  )
  for (file in progEstimFiles) {
    test_that(paste("plotBlup_1trait", basename(file)), {
      progBlup <- readProgBlupEstim(file)
      expect_no_error({
        trait <- names(progBlup[[1]][["blup_exp"]])
        p <- plotBlup_1trait(progBlup, trait = trait)
      })
      expect_true(all.equal(class(p), c("plotly", "htmlwidget")))
      expect_type(p$x$visdat[[1]]()$Cross, "character")

      for (s in c("asc", "dec", "???")) {
        expect_no_error({
          trait <- names(progBlup[[1]][["blup_exp"]])
          p <- plotBlup_1trait(progBlup, trait = trait)
        })
        expect_type(p$x$visdat[[1]]()$Cross, "character")
      }
    })
  }



  # plotBlup_2trait ----
  progEstimFiles <- c("../../data/results/progenyBlupEstimation_2traits.json")
  for (file in progEstimFiles) {
    test_that(paste("plotBlup_2trait", basename(file)), {
      progBlup <- readProgBlupEstim(file)
      expect_no_error({
        p <- plotBlup_2traits(
          blupDta = progBlup,
          x_trait = "trait1",
          y_trait = "trait2",
          confidenceLevel = 0.95,
          x_suffix = "",
          y_suffix = "",
          ellipses_npoints = 100
        )
      })
      expect_true(all.equal(class(p), c("plotly", "htmlwidget")))
    })
  }




  # evaluation_plot ----
  pheno_train <- readPhenoData("../../data/genomic_selection/pheno_train.csv")
  geno_train <- readGenoData("../../data/genomic_selection/geno_G1.vcf.gz")
  evaluation_results <- cross_validation_evaluation(
    pheno_train,
    geno_train,
    with_dominance = FALSE,
    n_folds = 3,
    n_repetitions = 2
  )
  test_that("evaluation_plot", {
    expect_no_error({
      p <- evaluation_plot(evaluation_results)
    })
    expect_true(all.equal(class(p), c("plotly", "htmlwidget")))
  })
})
