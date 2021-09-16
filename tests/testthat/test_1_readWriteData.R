# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test file for read/write data functions


capture_output({

  genoFiles <- c("../../data/geno/testMarkerData01.vcf",
                 "../../data/geno/testMarkerData01.vcf.gz")

  for (file in genoFiles) {
    test_that(paste("readGenoData", tools::file_ext(file)), {
      if (!file.exists(file)) {
        skip(paste("File", file, "not found. Skip test for this file"))
      }
      expect_error({
        gDta <- readGenoData(file)
      }, NA)
      expect_true(class(gDta) == "bed.matrix")

      # error
      expect_error({
        dta <- readGenoData("doNotExist")
      }, "File do not exists")
    })

    test_that(paste("downloadGenoData", tools::file_ext(file)), {
      if (!file.exists(file)) {
        skip(paste("File", file, "not found. Skip test for this file"))
      }
      file <- normalizePath(file)
      file <- paste0("file://", file)
      expect_error({
        gDta <- downloadGenoData(file)
      }, NA)
      expect_true(class(gDta) == "bed.matrix")

      # # error
      # expect_error({
      #   dta <- downloadGenoData("doNotExist")
      # })
    })
  }


  test_that("readPhenoData", {
    files <- c("../../data/pheno/testPhenoData01.csv")

    for (file in files) {
      expect_error({
        pDta <- readPhenoData(file)
      }, NA)
      expect_true(class(pDta) == "data.frame")
    }

    # error
    expect_error({
      dta <- readPhenoData("doNotExist")
    }, "File do not exists")
  })


  test_that("downloadPhenoData", {
    files <- c("../../data/pheno/testPhenoData01.csv")
    files <- normalizePath(files)
    files <- paste0("file://", files)
    for (file in files) {
      expect_error({
        pDta <- downloadPhenoData(file)
      }, NA)
      expect_true(class(pDta) == "data.frame")
    }

    # # error
    # expect_error({
    #   dta <- downloadPhenoData("doNotExist")
    # })
  })


  test_that("prepareData", {
    genoFiles <- c("../../data/geno/testMarkerData01.vcf.gz")
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
      }
    }

    # error
    expect_error({
      dta <- readData("doNotExist", phenfile)
    }, "File do not exists")
    expect_error({
      dta <- readData(genfile, "doNotExist")
    }, "File do not exists")
  })


  test_that("readData", {
    genoFiles <- c("../../data/geno/testMarkerData01.vcf.gz")
    phenoFiles <- c("../../data/pheno/testPhenoData01.csv")

    for (phenfile in phenoFiles) {
      for (genfile in genoFiles) {
        expect_error({
          dta <- readData(genfile, phenfile)
        }, NA)
        expect_true(class(dta$genoData) == "bed.matrix")
        expect_true(class(dta$phenoData) == "data.frame")
      }
    }

    # error
    expect_error({
      dta <- readData("doNotExist", phenfile)
    }, "File do not exists")
    expect_error({
      dta <- readData(genfile, "doNotExist")
    }, "File do not exists")
  })


  test_that("downloadData", {
    genoFiles <- c("../../data/geno/testMarkerData01.vcf.gz")
    phenoFiles <- c("../../data/pheno/testPhenoData01.csv")

    genoFiles <- normalizePath(genoFiles)
    genoFiles <- paste0("file://", genoFiles)

    phenoFiles <- normalizePath(phenoFiles)
    phenoFiles <- paste0("file://", phenoFiles)

    for (phenfile in phenoFiles) {
      for (genfile in genoFiles) {
        expect_error({
          dta <- downloadData(genfile, phenfile)
        }, NA)
        expect_true(class(dta$genoData) == "bed.matrix")
        expect_true(class(dta$phenoData) == "data.frame")
      }
    }

    # # error
    # expect_error({
    #   dta <- downloadData("doNotExist", phenfile)
    # })
    # expect_error({
    #   dta <- downloadData(genfile, "doNotExist")
    # })
  })


  test_that("readGWAS", {
    files <- c("../../data/results/gwasResult.json")

    for (file in files) {
      expect_error({
        gwas <- readGWAS(file)
      }, NA)
      expect_true(class(gwas) == "list")
      expect_true(all.equal(names(gwas),  c("gwas", "metadata")))
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
    }

    # error
    expect_error({
      dta <- readGWAS("doNotExist")
    }, "File do not exists")
  })


  test_that("downloadGWAS", {
    files <- c("../../data/results/gwasResult.json")
    files <- normalizePath(files)
    files <- paste0("file://", files)
    for (file in files) {
      expect_error({
        gwas <- downloadGWAS(file)
      }, NA)
      expect_true(class(gwas) == "list")
      expect_true(all.equal(names(gwas), c("gwas", "metadata")))
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
    }

    # # error
    # expect_error({
    #   dta <- downloadGWAS("doNotExist")
    # })
  })


  test_that("saveGWAS", {
    file <- c(g = "../../data/geno/testMarkerData01.vcf.gz",
              p = "../../data/pheno/testPhenoData01.csv")
    dta <- readData(file["g"], file["p"])
    resGwas <- gwas(dta,
                    trait = "Flowering.time.at.Arkansas",
                    test = "score",
                    fixed = 0,
                    response = "quantitative",
                    thresh_maf = 0.05,
                    thresh_callrate = 0.95)
    # this function is also tested in tests/testthat/test_2_gwas.R
    resfile <- tempfile()
    outfile <- saveGWAS(gwasRes = resGwas, metadata = "test", file = resfile)
    expect_equal(outfile, resfile)
    gwas <- readGWAS(outfile)
    expect_true(class(gwas) == "list")
    expect_true(all.equal(names(gwas), c("gwas", "metadata")))
    expect_true(class(gwas$gwas) == "data.frame")

    # error
    expect_error({
      file <- saveGWAS(gwasRes = resGwas, dir = "doNotExists")
    }, 'Error: "dir" directory should exists')
    expect_error({
      file <- saveGWAS(gwasRes = resGwas,  file = c(tempfile(), tempfile()))
    }, 'Error: only one file name should be provided')
  })
})
