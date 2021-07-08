# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test for gwas analysis

capture.output({
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
        if (length(na.omit(unique(dta$phenoData[,trait]))) > 2) {
          resp <- "quantitative"
        } else {
          if (test != "score") {
            next
          }
          resp <- "binary"
        }
        if (test %in% c("wald", "lrt")) {
          fix <- c(0:3)
        } else {fix <- 0}
        for (fixed in fix) {
          testName <- paste("gwas",resp,test,fixed,sep = "-")
          test_that(testName, {
            expect_error({
              resGwas <- gwas(dta,
                              trait,
                              test,
                              fixed = fixed,
                              response = resp,
                              thresh_maf = 0.05,
                              thresh_callrate = 0.95
              )
            }, NA)
            expect_true(class(resGwas) == "data.frame")
            expect_equal(colnames(resGwas),
                         resCols[[test]])
            expect_error({
              file <- saveGWAS(gwas = resGwas)
            }, "\"metadata\" is missing")
            expect_error({
              file <- saveGWAS(gwas = resGwas, metadata = list())
            }, NA)
            expect_error({
              gwas <- readGWAS(file)
            }, NA)
            expect_true(class(gwas$gwas) == "data.frame")
            expect_true(all.equal(gwas$gwas, resGwas))
            file.remove(file)
            expect_error({
              file <- saveGWAS(gwas = resGwas, dir = "doNotExists")
            }, 'Error: "dir" directory should exists')
          })
        }
      }
    }
  }



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
})
