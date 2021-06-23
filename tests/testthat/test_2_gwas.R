# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test for gwas analysis



test_that("gwas", {
  files <- list(c(g = "../../data/markers/testMarkerData01.vcf",
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
        } else { fix <- 0}
        for (fixed in fix) {
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
        }
      }
    }
  }
})
