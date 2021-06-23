# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test file for read/write data functions



test_that("readGenoData", {
  files <- c("../../data/markers/testMarkerData01.vcf",
             "../../data/markers/testMarkerData01.vcf.gz")

  for (file in files) {
    expect_error({
      gDta <- readGenoData(file)
    }, NA)
    expect_true(class(gDta) == "bed.matrix")
  }
})


test_that("readPhenoData", {
  files <- c("../../data/pheno/testPhenoData01.csv")

  for (file in files) {
    expect_error({
      pDta <- readPhenoData(file)
    }, NA)
    expect_true(class(pDta) == "data.frame")
  }
})


test_that("prepareData", {
  genoFiles <- c("../../data/markers/testMarkerData01.vcf")
  phenoFiles <- c("../../data/pheno/testPhenoData01.csv")

  for (phenfile in phenoFiles) {
    for (genfile in genoFiles) {
      gDta <- readGenoData(genfile)
      pDta <- readPhenoData(phenfile)
      expect_error({
        dta <- prepareData(gDta, pDta)
      }, NA)
      expect_true(class(dta$genoData) == "bed.matrix")
      expect_true(class(dta$phenoData) == "data.frame")
      expect_true(is.matrix(dta$grMatrix))
      expect_true(nrow(dta$phenoData) == nrow(dta$grMatrix))
    }
  }
})


test_that("readData", {
  genoFiles <- c("../../data/markers/testMarkerData01.vcf")
  phenoFiles <- c("../../data/pheno/testPhenoData01.csv")

  for (phenfile in phenoFiles) {
    for (genfile in genoFiles) {
      expect_error({
        dta <- readData(genfile,phenfile)
      }, NA)
      expect_true(class(dta$genoData) == "bed.matrix")
      expect_true(class(dta$phenoData) == "data.frame")
      expect_true(is.matrix(dta$grMatrix))
      expect_true(nrow(dta$phenoData) == nrow(dta$grMatrix))
    }
  }
})
