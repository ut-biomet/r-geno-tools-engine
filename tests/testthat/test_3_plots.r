# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Unit test for LD plots

capture_output({
  test_that("LD plot", {
    gDta <- readGenoData("../../data/markers/testMarkerData01.vcf")
    expect_error({
      imgFile <- LDplot(gDta, 1, 20, write = FALSE)
    }, NA)
    expect_error({
      imgFile <- LDplot(gDta, 42, 72, write = TRUE)
    }, NA)
    expect_error({
      imgFile <- LDplot(gDta, 1024, 1074, write = TRUE, dir = "testOutput")
    }, NA)
    file.remove(imgFile)

    expect_error({
      imgFile <- LDplot(gDta, 20, 2, write = FALSE)
    }, '"from" should be lower than "to"')
    expect_error({
      imgFile <- LDplot(gDta, 20, 200, write = FALSE)
    }, 'range size \\("to" - "from"\\) should be lower or equal than 50')
    expect_error({
      imgFile <- LDplot(gDta, 1024, 1074, write = TRUE, dir = "doNotExists")
    }, 'Error: "dir" directory should exists')
  })


  files <- list(c(g = "../../data/markers/testMarkerData01.vcf",
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
})
