# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Unit test for LD plots

capture_output({

  # LDplot
  test_that("LD plot", {
    gDta <- readGenoData("../../data/geno/testMarkerData01.vcf.gz")
    expect_error({
      imgFile <- LDplot(gDta, 1, 20, dir = NULL)
    }, NA)
    expect_error({
      imgFile <- LDplot(gDta, 42, 72, dir = tempdir())
    }, NA)
    expect_error({
      imgFile <- LDplot(gDta, 1024, 1074, dir = "testOutput")
    }, NA)
    file.remove(imgFile)

    # errors
    expect_error({
      imgFile <- LDplot(gDta, 1.5, 20, dir = NULL)
    },'`from` should be an positive integer of length 1' )
    expect_error({
      imgFile <- LDplot(gDta, 1, 20.1, dir = NULL)
    },'`to` should be an positive integer of length 1' )
    expect_error({
      imgFile <- LDplot(gDta, "2", 20, dir = NULL)
    }, '`from` should be an positive integer of length 1')
    expect_error({
      imgFile <- LDplot(gDta, 2, "20", dir = NULL)
    }, '`to` should be an positive integer of length 1')

    expect_error({
      imgFile <- LDplot(gDta, 20, 2, dir = NULL)
    }, '"from" should be lower than "to"')
    expect_error({
      imgFile <- LDplot(gDta, 20, 200, dir = NULL)
    }, 'range size \\("to" - "from"\\) should be lower or equal than 50')
    expect_error({
      imgFile <- LDplot(gDta, 1024, 1074, dir = "doNotExists")
    }, 'Error: "dir" directory should exists')
  })

  # manPlot
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
        for (chr in c(NA, unique(resGwas$chr))) {
          testName <- paste("manPlot", resp, test, paste0("chr:", chr),  sep = "-")
          test_that(testName,{
            expect_error({
              mP <- manPlot(gwas = resGwas,
                            adj_method = "bonferroni",
                            thresh_p = 0.5,
                            chr = chr,
                            title = "test manPlot")
            }, NA)
            expect_true(all.equal(class(mP), c("plotly", "htmlwidget")))
          })
        }
      }
    }
  }

  # Test manPlot() with wrong parameters:
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
                     chr = unique(resGwas$chr)[1])

  wrongParamsL <- list(p = list(c(runif(100, 0, 1), -0.1, 0.00015)),
                       adj_method = "do not exists",
                       thresh_p = list(-1, 1.02),
                       chr = c("doNotExists", 29))

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("manPlot WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          manPlot(gwas = params[["gwas"]],
                  adj_method = params[["adj_method"]],
                  thresh_p = params[["thresh_p"]],
                  chr = params[["chr"]])
        })
      })
    }
  }

})
