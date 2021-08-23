# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test for gwas analysis

capture.output({

  # Test gwas() with normal parameters
  files <- list(c(g = "../../data/geno/testMarkerData01.vcf.gz",
                  p = "../../data/pheno/testPhenoData01.csv"))
  resCols <- list("score" = c("chr", "pos", "id", "A1", "A2",
                              "freqA2", "score", "p"),
                  "wald" = c("chr", "pos", "id", "A1", "A2",
                             "freqA2", "h2", "beta", "sd","p"),
                  "lrt" = c("chr", "pos", "id", "A1", "A2",
                            "freqA2", "h2", "LRT", "p"))
  for (file in files) {
    dta <- readData(file["g"], file["p"])
    for (trait in c("Flowering.time.at.Arkansas", "Awn.presence")) {
      for (test in c("score", "wald", "lrt")) {
        if (length(na.omit(unique(dta$phenoData[, trait]))) > 2) {
          resp <- "quantitative"
        } else {
          resp <- "binary"
          if (test != "score") {
            next
          }
        }
        if (test %in% c("wald", "lrt")) {
          fix <- c(0:3)
        } else {fix <- 0}
        for (fixed in fix) {
          for (tMaf in seq(0,0.5,0.25)) {
            for (tCall in seq(0,1,0.5)) {
              testName <- paste("gwas",resp,test,fixed,tMaf,tCall,sep = "-")
              test_that(testName, {
                if (tMaf >= 0.5 || tCall >= 1) {
                  expect_warning({
                    resGwas <- gwas(dta,
                                    trait,
                                    test,
                                    fixed = fixed,
                                    response = resp,
                                    thresh_maf = tMaf,
                                    thresh_callrate = tCall)
                  })
                } else {
                  expect_error({
                    resGwas <- gwas(dta,
                                    trait,
                                    test,
                                    fixed = fixed,
                                    response = resp,
                                    thresh_maf = tMaf,
                                    thresh_callrate = tCall)
                  }, NA)
                }
                expect_true(class(resGwas) == "data.frame")
                expect_equal(colnames(resGwas), resCols[[test]])
                expect_error({
                  file <- saveGWAS(gwasRes = resGwas)
                }, "\"metadata\" is missing")
                expect_error({
                  file <- saveGWAS(gwasRes = resGwas, metadata = list())
                }, NA)
                expect_error({
                  resGwas2 <- readGWAS(file)
                }, NA)
                expect_true(class(resGwas2$gwas) == "data.frame")
                expect_true(all.equal(resGwas2$gwas, resGwas))
                file.remove(file)
              })
            }
          }
        }
      }
    }
  }


  # Test gwas() with wrong parameters:
  goodParams <- list(data = readData(files[[1]]["g"], files[[1]]["p"]),
                     trait = "Flowering.time.at.Arkansas",
                     test = "wald",
                     fixed = 0,
                     response = "quantitative",
                     thresh_maf = 0.05,
                     thresh_callrate = 0.95)

  wrongParamsL <- list(data = list(data.frame(c(1,2,3), c("a", "b", "c"))),
                       trait = list("Trait.that.do.not.exist", c("a","b")),
                       test = list("not.a.test", c("c", "d")),
                       fixed = list(-1, 1.5, "1", c(0,1)),
                       response = list("wrong resp", c("e", "f")),
                       thresh_maf = list(-0.06, 0.55, "0.1", c(0, 0.1)),
                       thresh_callrate = list(-0.02, 1.42, "0.21", c(0.85, 0.2)))

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("gwas, WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          gwas(data = params[["data"]],
               trait = params[["trait"]],
               test = params[["test"]],
               fixed = params[["fixed"]],
               response = params[["response"]],
               thresh_maf = params[["thresh_maf"]],
               thresh_callrate = params[["thresh_callrate"]])
        })
      })
    }
  }


  # Test adjustPval() with normal parameters
  s <- list(c(g = "../../data/geno/testMarkerData01.vcf.gz",
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
        for (adjMeth in c("holm",
                          "hochberg",
                          "bonferroni",
                          "BH",
                          "BY",
                          "fdr",
                          "none")) {
          for (thresh_p in c(NULL, 0.1, 0.05, 0.01)) {

            testName <- paste("adjustPval", resp, test, adjMeth, thresh_p, sep = "-")
            test_that(testName,{
              expect_error({
                adj <- adjustPval(p = resGwas$p,
                                  adj_method = adjMeth,
                                  thresh_p = thresh_p)
              }, NA)
              expect_true(is.list(adj))
              expect_equal(names(adj), c("p_adj", "thresh_adj"))
              if (adjMeth == "none") {
                expect_equal(thresh_p, adj$thresh_adj)
              }
            })
          }
        }
      }
    }
  }



  # Test adjustPval() with wrong parameters:
  goodParams <- list(p = runif(100, 0, 1),
                     adj_method = "bonferroni",
                     thresh_p = 0.05)

  wrongParamsL <- list(p = list(c(runif(100, 0, 1), -0.1, 0.00015)),
                       adj_method = "do not exists",
                       thresh_p = list(-1, 1.02))

  for (p in names(goodParams)) {
    wrongParams <- wrongParamsL[[p]]
    i <- 0
    for (wrongP in wrongParams) {
      i <- i+1
      params <- goodParams
      params[[p]] <- wrongP
      testName <- paste("adjustPval WrongParams", p, i, sep = "-")
      test_that(testName, {
        expect_error({
          adjustPval(p = params[["p"]],
                     adj_method = params[["adj_method"]],
                     thresh_p = params[["thresh_p"]])
        })
      })
    }
  }
})

