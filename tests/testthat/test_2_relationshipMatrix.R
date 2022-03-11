# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# relationship matrix tests.

capture_output({

  # pedRelMat -----
  pedFiles <- c('../../data/pedigree/testPedData_char.csv',
                '../../data/pedigree/testPedData_num.csv',
                '../../data/pedigree/testPedData_missFounder.csv')
  for (file in pedFiles) {
    test_that(paste("pedRelMat", basename(file)), {
      suppressWarnings({
        ped <- readPedData(file)
      })
      expect_error({
        relMat <- pedRelMat(ped)
      }, NA)
      expect_true(is.matrix(relMat))
      expect_true(is.numeric(relMat))
      expect_true(isSymmetric(relMat))
      expect_true(!any(is.na(relMat)))
      expect_true(all(relMat >= 0))
      expect_true(all(relMat <= 2))
      expect_true(!is.null(colnames(relMat)))
      expect_true(!is.null(row.names(relMat)))
      expect_identical(colnames(relMat), row.names(relMat))
      expect_identical(ped$data$ind, colnames(relMat))
    })
    test_that(paste("pedRelMat vs AGHmatrix", basename(file)),{
      skip_if_not_installed("AGHmatrix")
      suppressWarnings({
        ped <- readPedData(file)
      })
      exp_relMat <- AGHmatrix::Amatrix(replace(ped$data, is.na(ped$data), 0))
      relMat <- pedRelMat(ped)
      expect_identical(colnames(relMat), colnames(exp_relMat))
      expect_identical(row.names(relMat), row.names(exp_relMat))
      expect_identical(relMat, exp_relMat)
    })
  }


  # genoRelMat ----
  genoFiles <- c('../../data/geno/breedGame_geno.vcf.gz')

  for (file in genoFiles) {
    test_that(paste("genoRelMat", basename(file)), {
      suppressWarnings({
        geno <- readGenoData(file)
      })
      expect_error({
        relMat <- genoRelMat(geno)
      }, NA)
      expect_true(is.matrix(relMat))
      expect_true(is.numeric(relMat))
      expect_true(isSymmetric(relMat))
      expect_true(!any(is.na(relMat)))
      # expect_true(all(relMat >= 0)) # not for geno
      # expect_true(all(relMat <= 2)) # not for geno
      expect_true(!is.null(colnames(relMat)))
      expect_true(!is.null(row.names(relMat)))
      expect_identical(colnames(relMat), row.names(relMat))
      expect_identical(geno@ped$id, colnames(relMat))
    })
  }


})
