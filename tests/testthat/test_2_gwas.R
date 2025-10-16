# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test for gwas analysis

capture.output({
  # Test gwas() with normal parameters ----
  files <- list(c(
    g = "../../data/geno/testMarkerData01.vcf.gz",
    p = "../../data/pheno/testPhenoData01.csv"
  ))
  resCols <- list(
    "score" = c(
      "chr", "pos", "id", "A1", "A2",
      "freqA2", "score", "p"
    ),
    "wald" = c(
      "chr", "pos", "id", "A1", "A2",
      "freqA2", "h2", "beta", "sd", "p"
    ),
    "lrt" = c(
      "chr", "pos", "id", "A1", "A2",
      "freqA2", "h2", "LRT", "p"
    )
  )
  for (file in files) {
    for (tMaf in seq(0, 0.5, 0.25)) {
      for (tCall in seq(0, 1, 0.5)) {
        if (tCall == 1 & tMaf == 0.5) {
          # no markers matching those filters
          next
        }
        dta <- readData(file["g"], file["p"],
          maf_min = tMaf,
          callrate_min = tCall,
        )

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
            } else {
              fix <- 0
            }
            for (fixed in fix) {
              testName <- paste("gwas", resp, test, fixed, tMaf, tCall, sep = "-")
              test_that(testName, {
                expect_no_error({
                  resGwas <- gwas(dta,
                    trait,
                    test,
                    fixed = fixed,
                    response = resp
                  )
                })
                expect_true(class(resGwas) == "data.frame")
                expect_equal(colnames(resGwas), resCols[[test]])
                expect_no_error({
                  file <- saveGWAS(gwasRes = resGwas, metadata = list())
                })
                expect_no_error({
                  resGwas2 <- readGWAS(file)
                })
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

  # Test adjustPval() with normal parameters ----
  s <- list(c(
    g = "../../data/geno/testMarkerData01.vcf.gz",
    p = "../../data/pheno/testPhenoData01.csv"
  ))
  for (file in files) {
    dta <- readData(file["g"], file["p"],
      maf_min = 0.05,
      callrate_min = 0.95
    )
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
          response = resp
        )
        for (adjMeth in c(
          "holm",
          "hochberg",
          "bonferroni",
          "BH",
          "BY",
          "fdr",
          "none"
        )) {
          for (thresh_p in c(NULL, 0.1, 0.05, 0.01)) {
            for (empty in c(TRUE, FALSE)) {
              testName <- paste("adjustPval", resp, test, adjMeth, thresh_p,
                empty,
                sep = "-"
              )
              test_that(testName, {
                expect_no_error({
                  if (empty) {
                    p <- numeric()
                  } else {
                    p <- resGwas$p
                  }
                  adj <- adjustPval(
                    p = p,
                    adj_method = adjMeth,
                    thresh_p = thresh_p
                  )
                })
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
  }

  # adjustPval with edge case
  test_that("adjustPval handle gwas with some missing p values", {
    n_p <- 100
    p <- runif(n_p)

    missing_pval_indices <- sample.int(n_p, 5)
    p[missing_pval_indices] <- NA
    thresh_p <- 0.05

    for (adjMeth in c(
      "holm",
      "hochberg",
      "bonferroni",
      "BH",
      "BY",
      "fdr",
      "none"
    )) {
      expect_no_error({
        adj <- adjustPval(
          p = p,
          adj_method = adjMeth,
          thresh_p = thresh_p
        )
      })
      expect_true(all(is.na(adj$p_adj[missing_pval_indices])))
    }
  })

  test_that("adjustPval handle gwas with all missing p values", {
    p <- rep(NA, 100)
    thresh_p <- 0.05

    for (adjMeth in c(
      "holm",
      "hochberg",
      "bonferroni",
      "BH",
      "BY",
      "fdr",
      "none"
    )) {
      expect_no_error({
        adj <- adjustPval(
          p = p,
          adj_method = adjMeth,
          thresh_p = thresh_p
        )
      })
      expect_true(all(is.na(adj$p_adj)))
    }
  })
})
