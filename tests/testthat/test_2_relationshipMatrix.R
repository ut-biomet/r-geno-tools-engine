# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# relationship matrix tests.
library(testthat)

expect_correct_relmat <- function(relationship_matrix, individuals_names) {
  expect_true(is.matrix(relationship_matrix))
  expect_true(is.numeric(relationship_matrix))
  expect_true(isSymmetric(relationship_matrix))
  expect_true(!any(is.na(relationship_matrix)))
  expect_true(!is.null(colnames(relationship_matrix)))
  expect_true(!is.null(row.names(relationship_matrix)))
  expect_identical(
    colnames(relationship_matrix),
    row.names(relationship_matrix)
  )
  expect_identical(colnames(relationship_matrix), individuals_names)
}

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
      expect_no_error({
        relMat <- pedRelMat(ped)
      })
      expect_correct_relmat(relMat, ped$data$ind)
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
    test_that(paste("calc_additive_rel_mat", basename(file)), {
      geno <- readGenoData(file)
      # reduce `geno` size to be faster
      geno <- gaston::select.inds(geno, id %in% sample(geno@ped$id, 100))

      expect_no_error({
        relMat <- calc_additive_rel_mat(geno = geno, standardized = TRUE)
      })

      expect_true(is.list(relMat))
      expect_identical(names(relMat), c("geno_mat", "rel_mat"))

      expect_correct_relmat(relMat$rel_mat, geno@ped$id)

      expect_true(is.matrix(relMat$geno_mat))
      expect_true(is.numeric(relMat$geno_mat))
      expect_failure(expect_identical(relMat$geno_mat, gaston::as.matrix(geno)))
      expect_identical(
        calc_additive_rel_mat(geno = geno, standardized = FALSE)$geno_mat,
        gaston::as.matrix(geno)
      )
    })


    test_that(paste("calc_dominance_rel_mat", basename(file)), {
      geno <- readGenoData(file)
      # reduce `geno` size to be faster
      geno <- gaston::select.inds(geno, id %in% sample(geno@ped$id, 100))

      expect_no_error({
        relMat <- calc_dominance_rel_mat(geno = geno, standardized = TRUE)
      })

      expect_true(is.list(relMat))
      expect_identical(names(relMat), c("geno_mat", "rel_mat"))

      expect_correct_relmat(relMat$rel_mat, geno@ped$id)

      expect_true(is.matrix(relMat$geno_mat))
      expect_true(is.numeric(relMat$geno_mat))
      expect_failure(expect_identical(relMat$geno_mat, gaston::as.matrix(geno)))
      expect_identical(
        calc_additive_rel_mat(geno = geno, standardized = FALSE)$geno_mat,
        gaston::as.matrix(geno)
      )
    })


    test_that(paste("calc_additive_geno", basename(file)), {
      geno <- readGenoData(file)
      # reduce `geno` size to be faster
      geno <- gaston::select.inds(geno, id %in% sample(geno@ped$id, 100))

      expect_no_error({
        add_geno_std <- calc_additive_geno(geno = geno, standardized = TRUE)
        add_geno_not_std <- calc_additive_geno(geno = geno, standardized = FALSE)
      })

      for (add_geno in list(add_geno_std, add_geno_not_std)) {
        expect_true(is.matrix(add_geno))
        expect_true(is.numeric(add_geno))
      }

      expect_identical(add_geno_not_std, gaston::as.matrix(geno))
      expect_failure(expect_identical(add_geno_std, gaston::as.matrix(geno)))
    })


    test_that(paste("calc_dominance_geno", basename(file)), {
      geno <- readGenoData(file)
      # reduce `geno` size to be faster
      geno <- gaston::select.inds(geno, id %in% sample(geno@ped$id, 100))

      expect_no_error({
        dom_geno_std <- calc_dominance_geno(geno = geno, standardized = TRUE)
        dom_geno_not_std <- calc_dominance_geno(geno = geno, standardized = FALSE)
      })

      for (dom_geno in list(dom_geno_std, dom_geno_not_std)) {
        expect_true(is.matrix(dom_geno))
        expect_true(is.numeric(dom_geno))
      }
      expect_failure(expect_identical(dom_geno_std, dom_geno_not_std))

    })
  }


  test_that(paste("calc_dominance_geno missing values", basename(file)), {
      n = 1000
      p = 3
      geno <- matrix(sample(c(0, 1, 2), size = n*p, replace = TRUE),
                     ncol = p)
      geno[1, 1] <- NA
      geno <- gaston::as.bed.matrix(geno)
      dom_geno_std <- calc_dominance_geno(geno = geno, standardized = TRUE)
      expect_true(is.na(dom_geno_std[1, 1]))
      expect_true(!all(is.na(dom_geno_std[, 1])))
  })


  # combinedRelMat ----
  A <- pedRelMat(readPedData('../../data/pedigree/breedGame_pedigree.csv'))
  sampledInds <- sample(row.names(A), 100)
  A <- A[sampledInds, sampledInds]
  geno <- readGenoData('../../data/geno/breedGame_geno.vcf.gz')
  geno <- gaston::select.inds(geno, id %in% sampledInds)
  G <- calc_additive_rel_mat(geno = geno, standardized = TRUE)$rel_mat

  paramList <- list(
    default = list(ped_rm = A, geno_rm = G),
    legarra = list(ped_rm = A, geno_rm = G, method = 'Legarra'),
    martini = list(ped_rm = A, geno_rm = G, method = 'Martini',
                   tau = 0.5, omega = 0.5),
    martini_randParams = list(ped_rm = A, geno_rm = G, method = 'Martini',
                   tau = runif(1, 0, 20), omega = runif(1, -20, 1)),
    warn_missIndInPed = list(ped_rm = A[-(1:10), -(1:10)], geno_rm = G)
  )

  for (paramName in names(paramList)) {
    test_that(paste('combinedRelMat:', paramName), {

      params <- paramList[[paramName]]
      if (grepl('warn_', paramName)) {
        expect_warning({
          relMat <- do.call(combinedRelMat, params)
        })
      } else {
        expect_no_error({
          relMat <- do.call(combinedRelMat, params)
        })
      }
      expect_correct_relmat(relMat, colnames(params$ped_rm))
    })

    if (!grepl('warn_', paramName)) {
      test_that(paste("combineRelMat vs AGHmatrix:", paramName),{
        skip_if_not_installed("AGHmatrix")
        params <- paramList[[paramName]]
        suppressWarnings({
          relMat <- do.call(combinedRelMat, params)
        })
        agh_params <- list()
        agh_params$A <- params$ped_rm
        agh_params$G <- params$geno_rm
        agh_params$method <- "Martini"
        agh_params$omega <- params$omega
        agh_params$tau <- params$tau
        exp_relMat <- do.call(AGHmatrix::Hmatrix, agh_params)
        exp_relMat <- exp_relMat[row.names(agh_params$A),
                                 colnames(agh_params$A)]

        expect_identical(colnames(relMat), colnames(exp_relMat))
        expect_identical(row.names(relMat), row.names(exp_relMat))
        expect_true(all.equal.numeric(exp_relMat, relMat))
      })
    }
  }
})
