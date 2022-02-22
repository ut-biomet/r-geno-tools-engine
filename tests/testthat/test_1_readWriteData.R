# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# unit test file for read/write data functions


capture_output({

  genoFiles <- c("../../data/geno/testMarkerData01.vcf",
                 "../../data/geno/testMarkerData01.vcf.gz")

  for (file in genoFiles) {
    # readGenoData ----
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
      }, "Genotypic file do not exists")
    })

    # downloadGenoData ----
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


  # readPhenoData ----
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
    }, "Phenotypic file do not exists")

    expect_error({
      dta <- readPhenoData("../data/pheno_duplicated.csv")
    }, paste0("Duplicated individuals found in the phenotypic file. ",
              "Individuals must apprear only once in the phenotypic file."))
  })


  # downloadPhenoData ----
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

  # prepareData ----
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
    }, "Genotypic file do not exists")
    expect_error({
      dta <- readData(genfile, "doNotExist")
    }, "Phenotypic file do not exists")
  })


  # readData ----
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
    }, "Genotypic file do not exists")
    expect_error({
      dta <- readData(genfile, "doNotExist")
    }, "Phenotypic file do not exists")
  })


  # downloadData ----
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

  # readGWAS ----
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
    }, "GWAS file do not exists")
  })

  # downloadGWAS ----
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

  # saveGWAS ----
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







  # readPedData ----
  pedFiles <- c('../../data/pedigree/testPedData_char.csv',
                '../../data/pedigree/testPedData_num.csv',
                '../../data/pedigree/testPedData_missFounder.csv')
  for (file in pedFiles) {
    test_that(paste("readPedData", basename(file)), {
      if (grepl("missFounder", file)) {
        expect_warning({
          ped <- readPedData(file)
        })
      } else {
        expect_error({
          ped <- readPedData(file)
        }, NA)
      }
      expect_true(class(ped) == "list")
      expect_true(identical(names(ped), c("data", "graph")))
      expect_true(is.data.frame(ped$data))
      expect_true(identical(colnames(ped$data), c('ind', 'parent1', 'parent2')))
      expect_true(is.character(ped$data$ind))
      expect_true(is.character(ped$data$parent1))
      expect_true(is.character(ped$data$parent2))
      expect_true(identical(class(ped$graph), 'igraph'))

    })
  }
  # error
  test_that(paste("readPedData-error"), {
    expect_error({
      ped <- readPedData('doNotExist')
    }, "pedigree file do not exists")
    expect_error({
      ped <- readPedData('../data/pedigree_circleGenealogy.csv')
    }, 'inconsistent genealogy detected')
    expect_error({
      ped <- readPedData('../data/pedigree_duplicated.csv')
    }, 'inconsistent pedigree entrie\\(s\\)')
    expect_error({
      ped <- readPedData('../data/pedigree_NA.csv')
    }, 'unknown individual')
  })

  # downloadPedData ----
  for (file in pedFiles) {
    test_that(paste("downloadPedData",  basename(file)), {
      file <- normalizePath(file)
      file <- paste0("file://", file)
      if (grepl("missFounder", file)) {
        expect_warning({
          ped <- downloadPedData(file)
        })
      } else {
        expect_error({
          ped <- downloadPedData(file)
        }, NA)
      }

      expect_true(class(ped) == "list")
      expect_true(identical(names(ped), c("data", "graph")))
      expect_true(is.data.frame(ped$data))
      expect_true(identical(colnames(ped$data), c('ind', 'parent1', 'parent2')))
      expect_true(identical(class(ped$graph), 'igraph'))
    })
  }



  # saveRelMat ----
  for (file in pedFiles) {
    for (format in c('csv', 'json', 'notGood')) {
      test_that(paste("saveRelMat in", format, basename(file)), {
        if (format %in% c('csv', 'json')) {
          suppressWarnings({
            relMat <- pedRelMat(readPedData(file))
          })
          resfile <- tempfile()
          expect_error({
            outfile <- saveRelMat(relMat,
                                  metadata = list(info = 'test'),
                                  file = resfile,
                                  format = format)
          }, NA)
          expect_true(identical(outfile, resfile))
          relMat2 <- readRelMat(file = outfile, format = format)
          expect_true(identical(relMat2, relMat))

        } else {
          # error
          suppressWarnings({
            relMat <- pedRelMat(readPedData(file))
          })
          resfile <- tempfile()
          expect_error({
            outfile <- saveRelMat(relMat,
                                  metadata = list(info = 'test'),
                                  file = resfile,
                                  format = format)
          }, 'Error: File format misspecified. It should be either "csv",
               or "json".')
        }
      })
    }
  }





  # readRelMat ----
  relMatFiles <- c('../../data/results/pedigreeRelationship.csv',
                   '../../data/results/pedigreeRelationship.json')
  for (file in relMatFiles) {
    test_that(paste("readRelMat", basename(file)), {
      expect_error({
        relMat <- readRelMat(file)
      }, NA)
      expect_true(is.matrix(relMat))
      expect_true(isSymmetric(relMat))
      expect_true(!is.null(colnames(relMat)))
      expect_true(identical(colnames(relMat), row.names(relMat)))
    })
  }
  # error
  test_that(paste("readRelMat-error"), {
    expect_error({
      relMat <- readRelMat("doNotExist", format = 'csv')
    }, "Relationship matrix file do not exists")
  })


  # downloadRelMat ----
  for (file in relMatFiles) {
    test_that(paste("downloadRelMat",  basename(file)), {
      file <- normalizePath(file)
      file <- paste0("file://", file)
      expect_error({
        relMat <- downloadRelMat(file)
      }, NA)
      expect_true(is.matrix(relMat))
      expect_true(isSymmetric(relMat))
      expect_true(!is.null(colnames(relMat)))
      expect_true(identical(colnames(relMat), row.names(relMat)))
    })
  }

})
